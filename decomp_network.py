import networkx as nx

# Define decomposition network as a subclass of directed graph (DiGraph)
class decomp_network(nx.DiGraph):
    def __init__(self,pools=[],reactions=[]):
        # First initialize the network object
        super().__init__()
        for pool in pools:
            self.add_pool(**pool)
        for reaction in reactions:
            self.add_reaction(**reaction)
    def add_pool(self,name,initval=None,CN=None,immobile=True,**kwargs):
        self.add_node(name,immobile=immobile,**kwargs)
        if CN is not None:
            self.nodes[name]['CN']=CN
        elif 'CO2' not in name:
            print('Note: Pool %s has flexible CN ratio'%name)
        if immobile:
            if initval is not None:
                self.nodes[name]['initval']=initval
    def add_reaction(self,upstream_pool,downstream_pool,CUE,turnover_val,
                        turnover_units='y',turnover_name='TURNOVER_TIME',comment=None,
                        monod_terms=[],inhibition_terms=[],**kwargs):
        assert upstream_pool in self.nodes,'Pool %s must be added before reaction is added'%upstream_pool
        assert downstream_pool in self.nodes,'Pool %s must be added before reaction is added'%downstream_pool
        self.add_edge(upstream_pool,downstream_pool,CUE=CUE,turnover_val=turnover_val,turnover_units=turnover_units,**kwargs)
        if comment is not None:
            self.edges[(upstream_pool,downstream_pool)]['comment']=comment
        for term in monod_terms:
            for reqitem in ['species','k']:
                assert reqitem in term.keys(),'Monod terms must include "%s"'%reqitem
        self.edges[(upstream_pool,downstream_pool)]['monod_terms']=monod_terms
        for term in inhibition_terms:
            for reqitem in ['species','k','type']:
                assert reqitem in term.keys(),'Inhibition terms must include "%s"'%reqitem
        self.edges[(upstream_pool,downstream_pool)]['inhibition_terms']=inhibition_terms
        self.edges[(upstream_pool,downstream_pool)]['turnover_name']=turnover_name
            
        return
        
    def setup_writer(self):
        writer = PF_reaction_writer(self)
        return writer
    def write_Pflotran(self,**kwargs):
        return self.setup_writer().write_PF_reactions(**kwargs)
    def write_into_input_deck(self,templatefile_name,outputfile_name,**kwargs):
        return self.setup_writer().write_into_input_deck(templatefile_name,outputfile_name,**kwargs)
    def run_simulation(self,template_file,simulation_name,pflotran_exe,output_suffix='-obs-0.tec',print_output=False):
        inputdeck=simulation_name+'_generated.in'
        print('Setting up input deck in %s'%inputdeck)
        self.write_into_input_deck(template_file,inputdeck)
        import subprocess
        cmd='{pflotran_exe:s} -pflotranin {simname:s}_generated.in'.format(pflotran_exe=pflotran_exe,simname=simulation_name)
        print('Running cmd: %s'%cmd)
        result = subprocess.getoutput(cmd)
        if print_output:
            print(result)
        if 'ERROR' in result:
            print(result)
            raise RuntimeError('Pflotran simulation failed')
        import plot_pf_output
        outputfile=simulation_name + '_generated' + output_suffix
        # Set up for more flexibility in output formats
        output_data,units=plot_pf_output.read_tecfile(outputfile)
        return output_data,units
            


# Class for writing network out into Pflotran reaction sandbox format
class PF_reaction_writer:
    def __init__(self,network):
        self.network=network
        self.level=[]
        self.output=''
        self.indent_spaces=2
        return 

    def increase_level(self,levelname):
        self.add_line(levelname)
        self.level.append(levelname)
    def decrease_level(self):
        self.level.pop()
        self.output = self.output + ' '*(self.base_indent+len(self.level)*self.indent_spaces) + '/\n'
    def add_line(self,text):
        self.output = self.output + ' '*(self.base_indent+len(self.level)*self.indent_spaces) + text + '\n'
        
    def write_PF_reactions(self,base_indent=0,indent_spaces=2):
        # First write out all pools and their CN ratios
        self.indent_spaces=indent_spaces
        self.base_indent=base_indent

        # self.increase_level('REACTION_SANDBOX')
        self.increase_level('SOMDECOMP')
        self.increase_level('POOLS')
        for pool in self.network.nodes:
            if 'CO2' in pool:
                continue
            if 'CN' in self.network.nodes[pool]:
                self.add_line(pool.ljust(20)+'{CN:1.1f}'.format(CN=self.network.nodes[pool]['CN']))
            else:
                self.add_line(pool.ljust(20)+'# Variable C:N pool')
        self.decrease_level()
        for edge in self.network.edges:
            edge_data=self.network.edges[edge].copy()
            if 'comment' in edge_data.keys():
                self.add_line('# {comment:s}'.format(comment=edge_data.pop('comment')))
            self.increase_level('REACTION')
            self.add_line('UPSTREAM_POOL'.ljust(20)+edge[0])
            if 'CO2' not in edge[1]: # CO2 is included implicitly as the rest of CUE
                self.add_line('DOWNSTREAM_POOL'.ljust(20)+edge[1].ljust(20)+'{CUE:1.1e}'.format(CUE=edge_data.pop('CUE')))
            else:
                edge_data.pop('CUE',None)
            if 'turnover_name' in edge_data.keys():
                turnover_name=edge_data.pop('turnover_name')
            else:
                turnover_name='TURNOVER_TIME'
            self.add_line(turnover_name.ljust(20)+'{time:1.1e} {units:s}'.format(time=edge_data.pop('turnover_val'),units=edge_data.pop('turnover_units')))
            for monod in edge_data.pop('monod_terms'):
                    self.increase_level('MONOD')
                    self.add_line('SPECIES_NAME'.ljust(20)+monod['species'])
                    self.add_line('HALF_SATURATION_CONSTANT'+' {const:1.1e}'.format(const=monod['k']))
                    self.decrease_level()
            for inhib in edge_data.pop('inhibition_terms'):
                    self.increase_level('INHIBITION')
                    self.add_line('SPECIES_NAME'.ljust(20)+inhib['species'])
                    self.add_line('TYPE {inhtype:s}'.format(inhtype=inhib['type']))
                    self.add_line('INHIBITION_CONSTANT'+' {const:1.1e}'.format(const=inhib['k']))
                    self.decrease_level()
            # Write out the rest of the reaction attributes
            for param in edge_data.keys():
                val=edge_data[param]
                if isinstance(val,str):
                    self.add_line(param.ljust(20)+val)
                elif isinstance(val,float):
                    self.add_line(param.ljust(20)+' {val:1.2e}'.format(val=val))
                else:
                    raise ValueError('Parameter {param:s} was not a str or float'.format(param=param))
            self.decrease_level()
        
        for lev in range(len(self.level)):
            self.decrease_level()
        return self.output
        
    def write_into_input_deck(self,templatefile_name,outputfile_name,
            insert_reactions_text='[INSERT_SOMDECOMP_HERE]',
            insert_imspecies_text='[INSERT_IMMOBILE_SPECIES_HERE]',
            insert_constraints_text='[INSERT_IMMOBILE_CONSTRAINTS_HERE]',
            indent_spaces=2):
        base_indent=0
        with open(templatefile_name,'r') as templatefile:
            template_lines=templatefile.readlines()
        with open(outputfile_name,'w') as outputfile:
            for line in template_lines:
                if insert_reactions_text in line:
                    outputfile.write('\n' + ' '*(base_indent+indent_spaces) + '#### NOTE: Beginning of auto-inserted SOMDECOMP reactions####\n')
                    outputfile.write(self.write_PF_reactions(base_indent=base_indent+indent_spaces,indent_spaces=indent_spaces))
                    outputfile.write(' '*(base_indent+indent_spaces) + '#### NOTE: End of auto-inserted SOMDECOMP reactions ####\n\n')
                elif insert_imspecies_text in line:
                    outputfile.write('\n' + ' '*(base_indent+indent_spaces) + '#### NOTE: Beginning of auto-inserted immobile species ####\n')
                    for pool in self.network.nodes:
                        if not self.network.nodes[pool]['immobile']:
                            continue
                        if 'CN' in self.network.nodes[pool]:
                            outputfile.write(' '*(base_indent+indent_spaces) + pool + '\n')
                        else:
                            outputfile.write(' '*(base_indent+indent_spaces) + pool + 'C\n')
                            outputfile.write(' '*(base_indent+indent_spaces) + pool + 'N\n')
                    outputfile.write('\n' + ' '*(base_indent+indent_spaces) + '#### NOTE: End of auto-inserted immobile species ####\n')
                elif insert_constraints_text in line:
                    outputfile.write('\n' + ' '*(base_indent+indent_spaces) + '#### NOTE: Beginning of auto-inserted immobile species ####\n')
                    for pool in self.network.nodes:
                        if not self.network.nodes[pool]['immobile']:
                            continue
                        if 'CN' in self.network.nodes[pool]:
                            outputfile.write(' '*(base_indent+indent_spaces) + pool.ljust(20) + ' {const:1.1e}'.format(const=self.network.nodes[pool]['initval'])+'\n')
                        else:
                            if 'initCN' not in self.network.nodes[pool]:
                                raise ValueError('initCN must be provided for flexible CN pools')
                            outputfile.write(' '*(base_indent+indent_spaces) + (pool+'C').ljust(20) + '{const:1.1e}'.format(const=self.network.nodes[pool]['initval'])+'\n')
                            outputfile.write(' '*(base_indent+indent_spaces) + (pool+'N').ljust(20) + '{const:1.1e}'.format(const=self.network.nodes[pool]['initval']/self.network.nodes[pool]['initCN'])+'\n')
                    outputfile.write('\n' + ' '*(base_indent+indent_spaces) + '#### NOTE: End of auto-inserted immobile species ####\n')
                else:
                    if not line.isspace():
                        base_indent=len(line)-len(line.lstrip())
                    outputfile.write(line)
                    
        return


def make_nodecolors(nodes,POMcol='C0',microbecol='C1',DOMcol='C2',MAOMcol='C3',littercol='g',CWDcol='brown'):
    colors=[]
    for node in nodes:
        if 'MICROBES' in node:
            colors.append(microbecol)
        elif 'MAOM' in node or 'SOIL' in node:
            colors.append(MAOMcol)
        elif 'DOM' in node:
            colors.append(DOMcol)
        elif 'LITR' in node:
            colors.append(littercol)
        elif 'CWD' in node:
            colors.append(CWDcol)
        else:
            colors.append(POMcol)
    return colors
    
