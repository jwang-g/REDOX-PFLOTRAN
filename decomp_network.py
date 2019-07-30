import networkx as nx

class decomp_pool(dict):
    pass
    
class inhibition(dict):
    pass
    
class monod(dict):
    pass
    
class reaction(dict):
    def __init__(self,name,reactant_pools,product_pools,rate_constant,rate_units='y',reactiontype='SOMDECOMP',**kwargs):
        if not isinstance(reactant_pools,dict):
            raise TypeError('stoich must be a dictionary of pool names and stoichiometries')
        if not isinstance(product_pools,dict):
            raise TypeError('product_pools must be a dictionary of pool names and stoichiometries')
        self['name']=name
        self['reactant_pools']=reactant_pools
        self['product_pools']=product_pools
        self['rate_units']=rate_units
        if reactiontype not in ['MICROBIAL','SOMDECOMP']:
            raise ValueError('Only MICROBIAL and SOMDECOMP reactions are currently implemented')
        self['reactiontype']=reactiontype
        self['rate_constant']=rate_constant
        self.update(kwargs)
        


# Class for writing network out into Pflotran reaction sandbox format
class PF_writer:
    def __init__(self,network,indent_spaces=2,base_indent=0):
        self.level=[]
        self.output=''
        self.indent_spaces=indent_spaces
        self.base_indent=base_indent
        self.network=network
        return 
        
    def reset(self):
        self.__init__(self.network)

    def current_indent(self):
        return self.base_indent+len(self.level)
    def increase_level(self,levelname):
        self.add_line(levelname)
        self.level.append(levelname)
    def decrease_level(self):
        self.level.pop()
        self.output = self.output + ' '*(self.base_indent+len(self.level)*self.indent_spaces) + '/\n'
    def add_line(self,text):
        self.output = self.output + ' '*(self.base_indent+len(self.level)*self.indent_spaces) + text + '\n'
    def write(self,*args,**kwargs):
        raise NotImplementedError('Must use subclass of PF_reaction_writer to write out reactions')
        
class PF_network_writer(PF_writer):
    def write_all_reactions(self,base_indent=0,indent_spaces=2):
        self.indent_spaces=indent_spaces
        self.base_indent=base_indent
        
        # Microbial reactions
        already_done = []
        for (ins,outs,reaction) in self.network.edges(data='reaction'):
            if reaction['reactiontype'] == 'MICROBIAL':
                if reaction['name'] in already_done:
                    continue
                already_done.append(reaction['name'])
                self.output = self.output + PF_microbial_reaction_writer(reaction,base_indent=self.current_indent()).write_reaction()
        
        # List all sandbox reactions
        sandbox_reacts=self.network.edge_subgraph([(ins,outs,ks) for ins,outs,ks,r in self.network.edges(data='reaction',keys=True) if r['reactiontype']=='SOMDECOMP'])
        # First write out all pools and their CN ratios
        if len(sandbox_reacts)>0:
            self.increase_level('REACTION_SANDBOX')
            self.increase_level('SOMDECOMP')
            self.increase_level('POOLS')
            for pool in sandbox_reacts.nodes:
                if 'CO2' in pool or 'HCO3-' in pool:
                    continue
                if 'CN' in sandbox_reacts.nodes[pool]:
                    self.add_line(pool.ljust(20)+'{CN:1.1f}'.format(CN=sandbox_reacts.nodes[pool]['CN']))
                else:
                    self.add_line(pool.ljust(20)+'# Variable C:N pool')
            self.decrease_level()
        
            # Write out all the reactions
            already_done=[]
            CO2name = None
            for ins,outs,reaction in sandbox_reacts.edges(data='reaction'):
                if reaction['name'] in already_done:
                    continue
                already_done.append(reaction['name'])
                self.output = self.output + PF_sandbox_reaction_writer(reaction,base_indent=self.current_indent()).write_reaction()
                # Figure out which CO2 species is used in these reactions. Assumes all reactions use the same one so any can be picked
                for spec in ['HCO3-','CO2(g)','CO2(aq)','CO2(g)*']:
                    if spec in reaction['product_pools']:
                        CO2name=spec
                        
            if CO2name is None:
                print('WARNING: CO2 name not found in SOMDECOMP reactions')
            else:
                self.add_line('CO2_SPECIES_NAME '+CO2name)
            
        for lev in range(len(self.level)):
            self.decrease_level()
        return self.output
            


    def write_into_input_deck(self,templatefile_name,outputfile_name,constraintname='initial',
            indent_spaces=2,length_days=None):
        base_indent=0
        with open(templatefile_name,'r') as templatefile:
            template_lines=templatefile.readlines()
    
        for line in template_lines:
            if 'CHEMISTRY' in line and not line.strip().startswith('#'):
                self.output = self.output + line
                # Start of chemistry block
                # Primary species
                self.increase_level('PRIMARY_SPECIES')
                self.add_line('#### NOTE: Beginning of auto-inserted primary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='primary':
                        self.add_line( pool )
                self.add_line( '#### NOTE: End of auto-inserted primary species ####')
                self.decrease_level()
                
                # Assume that all primary species should be put into decoupled eq reactions just in case
                self.increase_level('DECOUPLED_EQUILIBRIUM_REACTIONS')
                self.add_line( '#### NOTE: Beginning of auto-inserted primary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='primary':
                        self.add_line( pool )
                self.add_line( '#### NOTE: End of auto-inserted primary species ####')
                self.decrease_level()
                
                # Secondary species
                self.increase_level('SECONDARY_SPECIES')
                self.add_line('#### NOTE: Beginning of auto-inserted secondary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='secondary':
                        self.add_line( pool )
                self.add_line( '#### NOTE: End of auto-inserted secondary species ####')
                self.decrease_level()
                
                # Minerals
                self.increase_level('MINERALS')
                self.add_line( '#### NOTE: Beginning of auto-inserted secondary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='mineral':
                        self.add_line( pool)
                self.add_line( '#### NOTE: End of auto-inserted secondary species ####')
                self.decrease_level()
                
                # Mineral kinetics
                self.increase_level('MINERAL_KINETICS')
                self.add_line( '#### NOTE: Beginning of auto-inserted secondary species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='mineral':
                        self.increase_level(pool)
                        self.add_line('RATE_CONSTANT  ' + self.network.nodes[pool]['rate'])
                        self.decrease_level()
                self.add_line( '#### NOTE: End of auto-inserted secondary species ####')
                self.decrease_level()
                
                # Immobile species
                self.increase_level('IMMOBILE_SPECIES')
                self.add_line( '#### NOTE: Beginning of auto-inserted immobile species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='immobile':
                        if 'CN' in self.network.nodes[pool] or pool in ['HRimm','Nmin','Nimp','Nimm','NGASmin']:
                            self.add_line( pool)
                        else:
                            self.add_line( pool + 'C')
                            self.add_line( pool + 'N')
                self.add_line( '#### NOTE: End of auto-inserted immobile species ####')
                self.decrease_level()
                
                # Gas species
                self.increase_level('GAS_SPECIES')
                self.add_line( '#### NOTE: Beginning of auto-inserted gas species ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='gas':
                        self.add_line( pool)
                self.add_line( '#### NOTE: End of auto-inserted gas species ####')
                self.decrease_level()
                
                # Reactions
                self.add_line( '#### NOTE: Beginning of auto-inserted reactions ####')
                self.write_all_reactions(base_indent=base_indent+indent_spaces,indent_spaces=indent_spaces)
                self.add_line( '#### NOTE: End of auto-inserted reactions ####')
            
            elif line.strip().startswith('CONSTRAINT') and line.strip().endswith(constraintname):
                self.output = self.output + line
                self.increase_level('IMMOBILE')
                self.add_line( '#### NOTE: Beginning of auto-inserted immobile species ####')
                for pool in self.network.nodes:
                    if not self.network.nodes[pool]['kind']=='immobile':
                        continue
                    if 'CN' in self.network.nodes[pool] or pool in ['HRimm','Nmin','Nimp','Nimm','NGASmin']:
                        self.add_line( pool.ljust(20) + ' {const:1.1e}'.format(const=self.network.nodes[pool]['initval']))
                    else:
                        if 'initCN' not in self.network.nodes[pool]:
                            raise ValueError('initCN must be provided for flexible CN pools: pool %s'%pool)
                        self.add_line( (pool+'C').ljust(20) + '{const:1.1e}'.format(const=self.network.nodes[pool]['initval']))
                        self.add_line( (pool+'N').ljust(20) + '{const:1.1e}'.format(const=self.network.nodes[pool]['initval']/self.network.nodes[pool]['initCN']))
                self.add_line( '#### NOTE: End of auto-inserted immobile species ####')
                self.decrease_level()
                
                self.increase_level('CONCENTRATIONS')
                self.add_line( '#### NOTE: Beginning of auto-inserted concentration constraints ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='primary':
                        if isinstance(self.network.nodes[pool]['initval'],str):
                            self.add_line( (pool).ljust(20) + self.network.nodes[pool]['initval'])
                        else:
                            self.add_line( (pool).ljust(20) + '{const:1.1e}'.format(const=self.network.nodes[pool]['initval']))
                            
                self.add_line( '#### NOTE: End of auto-inserted concentration constraints ####')
                self.decrease_level()
                
                self.increase_level('MINERALS')
                self.add_line( '#### NOTE: Beginning of auto-inserted mineral constraints ####')
                for pool in self.network.nodes:
                    if self.network.nodes[pool]['kind']=='mineral':
                        if isinstance(self.network.nodes[pool]['initval'],str):
                            self.add_line( (pool).ljust(20) + self.network.nodes[pool]['initval'])
                        else:
                            self.add_line( (pool).ljust(20) + '{const:1.1e}'.format(const=self.network.nodes[pool]['initval']))
                            
                self.add_line( '#### NOTE: End of auto-inserted mineral constraints ####')
                self.decrease_level()
                
            elif 'FINAL_TIME' in line and length_days is not None:
                self.add_line( 'FINAL_TIME {ndays:1.1e} d\n'.format(ndays=length_days))
            else:
                if not line.isspace():
                    base_indent=len(line)-len(line.lstrip())
                self.output = self.output + line
                    
        with open(outputfile_name,'w') as outputfile:
            outputfile.write(self.output)
        return

        
        
    def run_simulation(self,template_file,simulation_name,pflotran_exe,output_suffix='-obs-0.tec',print_output=False,length_days=None):
        inputdeck=simulation_name+'_generated.in'
        print('Setting up input deck in %s'%inputdeck)
        self.write_into_input_deck(template_file,inputdeck,length_days=length_days)
        import subprocess
        cmd='{pflotran_exe:s} -pflotranin {simname:s}_generated.in'.format(pflotran_exe=pflotran_exe,simname=simulation_name)
        print('Running cmd: %s'%cmd)
        status,output = subprocess.getstatusoutput(cmd)
        print('Simulation finished with status %d'%status)
        if print_output:
            print(output)
        if 'ERROR' in output or 'Stopping' in output:
            print(output)
            raise RuntimeError('Pflotran simulation failed')
        import plot_pf_output
        outputfile=simulation_name + '_generated' + output_suffix
        # Set up for more flexibility in output formats
        output_data,units=plot_pf_output.read_tecfile(outputfile)
        return output_data,units
        
class PF_sandbox_reaction_writer(PF_writer):
    def write_reaction(self,base_indent=None,indent_spaces=None):
        reaction_data=self.network.copy()
        assert len(reaction_data['reactant_pools']) == 1,'Only one reactant pool is allowed for sandbox reactions'
        assert list(reaction_data['reactant_pools'].values())[0] == 1.0,'Stoichiometry of upstream pool in sandbox reaction must be 1.0'
        if base_indent is not None:
            self.base_indent=base_indent
        self.add_line('# {name:s}'.format(name=reaction_data.pop('name')))
        if 'comment' in reaction_data.keys():
            self.add_line('# {comment:s}'.format(comment=reaction_data.pop('comment')))
        self.increase_level('REACTION')
        upstream_pool=list(reaction_data.pop('reactant_pools').keys())[0]
        self.add_line('UPSTREAM_POOL'.ljust(20)+upstream_pool)
        downstream_pools = reaction_data.pop('product_pools')
        for pool in downstream_pools.keys():
            if 'CO2' not in pool and 'HCO3-' not in pool: # CO2 is included implicitly as the rest of CUE
                self.add_line('DOWNSTREAM_POOL'.ljust(20)+pool.ljust(20)+'{CUE:1.1e}'.format(CUE=downstream_pools[pool]))
        turnover_name=reaction_data.pop('turnover_name','TURNOVER_TIME')
        self.add_line(turnover_name.ljust(20)+'{time:1.1e} {units:s}'.format(time=reaction_data.pop('rate_constant'),units=reaction_data.pop('rate_units')))
        for monod in reaction_data.pop('monod_terms',[]):
                self.increase_level('MONOD')
                self.add_line('SPECIES_NAME'.ljust(20)+monod['species'])
                self.add_line('HALF_SATURATION_CONSTANT'+' {const:1.1e}'.format(const=monod['k']))
                if 'threshold' in monod:
                    self.add_line('THRESHOLD_CONCENTRATION {const:1.1e}'.format(const=monod['threshold']))
                self.decrease_level()
        for inhib in reaction_data.pop('inhibition_terms',[]):
                self.increase_level('INHIBITION')
                self.add_line('SPECIES_NAME'.ljust(20)+inhib['species'])
                self.add_line('TYPE {inhtype:s}'.format(inhtype=inhib['type']))
                self.add_line('INHIBITION_CONSTANT'+' {const:1.1e}'.format(const=inhib['k']))
                self.decrease_level()
        # Write out the rest of the reaction attributes
        reaction_data.pop('reactiontype',None)
        for param in reaction_data.keys():
            val=reaction_data[param]
            if isinstance(val,str):
                self.add_line(param.ljust(20)+val)
            elif isinstance(val,float):
                self.add_line(param.ljust(20)+' {val:1.2e}'.format(val=val))
            else:
                raise ValueError('Parameter {param:s} was not a str or float'.format(param=param))
        self.decrease_level()
    
        return self.output
        
        
class PF_microbial_reaction_writer(PF_writer):
    def write_reaction(self,base_indent=None,indent_spaces=None):
        reaction_data=self.network.copy()
        
        if base_indent is not None:
            self.base_indent=base_indent
        self.add_line('# {name:s}'.format(name=reaction_data.pop('name')))
        if 'comment' in reaction_data.keys():
            self.add_line('# {comment:s}'.format(comment=reaction_data.pop('comment')))
        self.increase_level('MICROBIAL_REACTION')
        line = 'REACTION '
        reactants=reaction_data.pop('reactant_pools')
        products =reaction_data.pop('product_pools')
        first = True
        for chem in reactants:
            if first:
                first = False
            else:
                line = line + ' + '
            line = line + '{stoich:1.1e} {chem:s} '.format(stoich=reactants[chem],chem=chem)
        line = line + ' -> '
        first = True
        for chem in products:
            if first:
                first = False
            else:
                line = line + ' + '
            line = line + '{stoich:1.1e} {chem:s} '.format(stoich=products[chem],chem=chem)
        self.add_line(line)

        self.add_line('RATE_CONSTANT'.ljust(20)+'{time:1.1e}'.format(time=reaction_data.pop('rate_constant')))
        for monod in reaction_data.pop('monod_terms',[]):
                self.increase_level('MONOD')
                self.add_line('SPECIES_NAME'.ljust(20)+monod['species'])
                self.add_line('HALF_SATURATION_CONSTANT'+' {const:1.1e}'.format(const=monod['k']))
                if 'threshold' in monod:
                    self.add_line('THRESHOLD_CONCENTRATION {const:1.1e}'.format(const=monod['threshold']))
                self.decrease_level()
        for inhib in reaction_data.pop('inhibition_terms',[]):
                self.increase_level('INHIBITION')
                self.add_line('SPECIES_NAME'.ljust(20)+inhib['species'])
                self.add_line('TYPE {inhtype:s}'.format(inhtype=inhib['type']))
                self.add_line('INHIBITION_CONSTANT'+' {const:1.1e}'.format(const=inhib['k']))
                self.decrease_level()
        # Write out the rest of the reaction attributes
        # What to do with biomass? Maybe better not to use it
        self.decrease_level()
    
        return self.output
        

# Define decomposition network as a subclass of directed graph (DiGraph)
class decomp_network(nx.MultiDiGraph):
    def __init__(self,pools=[],reactions=[]):
        # First initialize the network object
        super().__init__()
        for pool in pools:
            self.add_pool(pool)
        for reaction in reactions:
            self.add_reaction(reaction)
    def add_pool(self,pool):
        pool_data=pool.copy()
        CN=pool_data.pop('CN',None)
        name=pool_data.pop('name')
        kind=pool_data.pop('kind')
        initval=pool_data.pop('initval',1e-15)
        self.add_node(name,kind=kind,**pool_data)
        if CN is not None:
            self.nodes[name]['CN']=CN
        self.nodes[name]['initval']=initval
    def add_reaction(self,reaction):
        reaction_data=reaction.copy()
        reactant_pools=reaction_data.pop('reactant_pools')
        product_pools=reaction_data.pop('product_pools')
        for upstream_pool in reactant_pools.keys():
            assert upstream_pool in self.nodes,'Pool %s must be added before reaction is added'%upstream_pool
            for downstream_pool in product_pools.keys():
                assert downstream_pool in self.nodes,'Pool %s must be added before reaction is added'%downstream_pool
                monod_terms=reaction_data.pop('monod_terms',[])
                inhibition_terms=reaction_data.pop('inhibition_terms',[])
                k=self.add_edge(upstream_pool,downstream_pool,name=reaction_data['name'],reaction=reaction)
                for term in monod_terms:
                    for reqitem in ['species','k']:
                        assert reqitem in term.keys(),'Monod terms must include "%s"'%reqitem
                self.edges[(upstream_pool,downstream_pool,k)]['monod_terms']=monod_terms
                for term in inhibition_terms:
                    for reqitem in ['species','k','type']:
                        assert reqitem in term.keys(),'Inhibition terms must include "%s"'%reqitem
                self.edges[(upstream_pool,downstream_pool,k)]['inhibition_terms']=inhibition_terms
            
        return
        


def make_nodecolors(nodes,POMcol='C0',microbecol='C1',DOMcol='C2',MAOMcol='C3',littercol='g',CWDcol='brown',
                            mineralcol='C4',primarycol='C5',secondarycol='C7',gascol='C7'):
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
        elif 'kind' in nodes[node]:
            if nodes[node]['kind']=='mineral':
                colors.append(mineralcol)
            elif nodes[node]['kind']=='primary':
                colors.append(primarycol)
            elif nodes[node]['kind']=='secondary':
                colors.append(secondarycol)
            elif nodes[node]['kind']=='gas':
                colors.append(gascol)
            else:
                 colors.append(primarycol)
        else:
            colors.append(POMcol)
    return colors
    
def get_reaction_from_database(filename,pool):
    with open(filename,'r') as dbase:
        for line in dbase:
            if line.startswith("'%s'"%pool['name']):
                out=nx.MultiDiGraph()
                lsplit=line.split()
                if pool['kind']=='secondary':
                    offset=1
                elif pool['kind']=='mineral':
                    offset=2
                elif pool['kind']=='gas':
                    offset=2
                else:
                    raise TypeError('Pool must be secondary, mineral, or gas')
                nspecies=int(lsplit[offset])
                for specnum in range(nspecies):
                    out.add_edge(pool['name'],lsplit[offset+2+specnum*2].strip("'"),reactiontype='equilibrium')
                    out.add_edge(lsplit[offset+2+specnum*2].strip("'"),pool['name'],reactiontype='equilibrium')
                return out
        raise ValueError('Species %s not found in database'%pool['name'])

