import networkx as nx

# Simple decomposition network including one POM pool, one microbial pool, one SOM pool derived from microbial necromass, and one DOM pool
decomp_network_simple=nx.DiGraph()

# Pools are nodes, with attributes that control their properties (e.g. C:N ratio)
# Could add more things like initial conditions if we want to generate a whole input deck
decomp_network_simple.add_node('POM',CN=15)
decomp_network_simple.add_node('SOIL_MICROBES',CN=10)
decomp_network_simple.add_node('SOM2',CN=10)
decomp_network_simple.add_node('DOM1',CN=10)

# Decomposition pathways are edges (connecting two nodes) with attributes that control the reaction
decomp_network_simple.add_edge('POM','SOIL_MICROBES',CUE=0.6,turnover_val=0.5,turnover_units='y',comment='Microbial decomposition of POM',
                                monod_species='SOIL_MICROBES',monod_k=0.1e-3)
decomp_network_simple.add_edge('SOIL_MICROBES','SOM2',CUE=0.1,turnover_val=0.1,turnover_units='y',comment='Microbial biomass turnover to SOM2',
                                inhibition_species='SOIL_MICROBES',inhibition_type='INVERSE_MONOD',inhibition_k=1e-11)
decomp_network_simple.add_edge('SOM2','SOIL_MICROBES',CUE=0.6,turnover_val=10,turnover_units='y',comment='Microbial decomposition of SOM2',
                                monod_species='SOIL_MICROBES',monod_k=0.1e-3)
decomp_network_simple.add_edge('POM','DOM1',CUE=1.0,turnover_val=0.1,turnover_units='y',comment='Dissolution of POM into DOM1',
                                inhibition_species='DOM1',inhibition_type='MONOD',inhibition_k=0.1e-12)


# A more complex decomposition network approximating the structure of the CORPSE model

decomp_network_CORPSE=nx.DiGraph()
decomp_network_CORPSE.add_node('LABILE_POM',CN=    15.0)
decomp_network_CORPSE.add_node('RESISTANT_POM',CN= 15.0)
decomp_network_CORPSE.add_node('DOM1',CN=    15.0)
decomp_network_CORPSE.add_node('DOM2',CN= 15.0)
decomp_network_CORPSE.add_node('SOIL_MICROBES',CN= 10.0)
decomp_network_CORPSE.add_node('NECROMASS',CN=     10.0)
decomp_network_CORPSE.add_node('DOM3',CN= 10.0)
decomp_network_CORPSE.add_node('MAOM',CN=          10.0)

# Microbial decomposition of unprotected SOM
decomp_network_CORPSE.add_edge('LABILE_POM','SOIL_MICROBES',CUE=0.6,turnover_val=0.1,turnover_units='y',comment='Labile POM microbial decomposition',
                                monod_species='SOIL_MICROBES',monod_k=0.1e-3)
decomp_network_CORPSE.add_edge('RESISTANT_POM','SOIL_MICROBES',CUE=0.1,turnover_val=0.5,turnover_units='y',comment='Resistant POM microbial decomposition',
                                monod_species='SOIL_MICROBES',monod_k=0.1e-3)
decomp_network_CORPSE.add_edge('NECROMASS','SOIL_MICROBES',CUE=0.6,turnover_val=0.1,turnover_units='y',comment='Necromass microbial decomposition',
                                monod_species='SOIL_MICROBES',monod_k=0.1e-3)
                                
# Microbial biomass turnover. Here, CUE would determine maintenance respiration loss as a fraction of turnover
decomp_network_CORPSE.add_edge('SOIL_MICROBES','NECROMASS',CUE=0.6,turnover_val=0.1,turnover_units='y',
                                comment='Microbial biomass turnover. Inhibition  prevents biomass from declining to zero',
                                inhibition_species='SOIL_MICROBES',inhibition_type='INVERSE_MONOD',inhibition_k=0.1e-11)
                                
# Physical/mineral protection of necromass. CUE=1.0 because it is an abiotic conservative reaction
decomp_network_CORPSE.add_edge('NECROMASS','MAOM',CUE=1.0,turnover_val=0.01,turnover_units='y',comment='Physical/mineral protection of necromass',
                                inhibition_species='MAOM',inhibition_k=1e-7,inhibition_type='MONOD')
# First order turnover of MAOM back to necromass
decomp_network_CORPSE.add_edge('MAOM','NECROMASS',CUE=1.0,turnover_val=50.0,turnover_units='y',comment='First order turnover of MAOM back to necromass') 

# Dissolution into DOM
decomp_network_CORPSE.add_edge('LABILE_POM','DOM1',CUE=1.0,turnover_val=0.1,turnover_units='y',comment='Dissolution of labile POM into DOM',
                                inhibition_species='DOM1',inhibition_type='MONOD',inhibition_k=1e-12)
decomp_network_CORPSE.add_edge('RESISTANT_POM','DOM2',CUE=1.0,turnover_val=0.1,turnover_units='y',comment='Dissolution of resistant POM into DOM',
                                inhibition_species='DOM2',inhibition_type='MONOD',inhibition_k=1e-12)
decomp_network_CORPSE.add_edge('NECROMASS','DOM3',CUE=1.0,turnover_val=0.1,turnover_units='y',comment='Dissolution of necromass into DOM',
                                inhibition_species='DOM3',inhibition_type='MONOD',inhibition_k=1e-12)      

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

        self.increase_level('REACTION_SANDBOX')
        self.increase_level('SOMDECOMP')
        self.increase_level('POOLS')
        for pool in self.network.nodes:
            self.add_line(pool.ljust(20)+'{CN:1.1f}'.format(CN=self.network.nodes[pool]['CN']))
        self.decrease_level()
        for edge in self.network.edges:
            edge_data=self.network.edges[edge]
            if 'comment' in edge_data.keys():
                self.add_line('# {comment:s}'.format(comment=edge_data['comment']))
            self.increase_level('REACTION')
            self.add_line('UPSTREAM_POOL'.ljust(20)+edge[0])
            self.add_line('DOWNSTREAM_POOL'.ljust(20)+edge[1].ljust(20)+'{CUE:1.1e}'.format(CUE=edge_data['CUE']))
            self.add_line('TURNOVER_TIME'.ljust(20)+'{time:1.1e} {units:s}'.format(time=edge_data['turnover_val'],units=edge_data['turnover_units']))
            if 'monod_species' in edge_data.keys():
                if edge_data['monod_species'] is None:
                    pass
                elif isinstance(edge_data['monod_species'],str):
                    self.increase_level('MONOD')
                    self.add_line('SPECIES_NAME'.ljust(20)+edge_data['monod_species'])
                    self.add_line('HALF_SATURATION_CONSTANT'+' {const:1.1e}'.format(const=edge_data['monod_k']))
                    self.decrease_level()
                else: # Assume it is a list
                    for num in range(len(edge_data['monod_species'])):
                        self.add_level('MONOD')
                        self.add_line('SPECIES_NAME'.ljust(20)+edge_data['monod_species'][num])
                        self.add_line('HALF_SATURATION_CONSTANT'+' {const:1.1e}'.format(const=edge_data['monod_k'][num]))
                        self.decrease_level()
            if 'inhibition_species' in edge_data.keys():
                if edge_data['inhibition_species'] is None:
                    pass
                elif isinstance(edge_data['inhibition_species'],str):
                    self.increase_level('INHIBITION')
                    self.add_line('SPECIES_NAME'.ljust(20)+edge_data['inhibition_species'])
                    self.add_line('TYPE {inhtype:s}'.format(inhtype=edge_data['inhibition_type']))
                    self.add_line('INHIBITION_CONSTANT'+' {const:1.1e}'.format(const=edge_data['inhibition_k']))
                    self.decrease_level()
                else: # Assume it is a list
                    for num in range(len(edge_data['inhibition_species'])):
                        self.increase_level('INHIBITION')
                        self.add_line('SPECIES_NAME'.ljust(20)+edge_data['inhibition_species'][num])
                        self.add_line('TYPE {inhtype:s}'.format(edge_data['inhibition_type'][num]))
                        self.add_line('INHIBITION_CONSTANT'+' {const:1.1e}'.format(edge_data['inhibition_k'][num]))
                        self.decrease_level()
            self.decrease_level()
        
        for lev in range(len(self.level)):
            self.decrease_level()
        return self.output

print('\n\n#### Simple network\n\n')
print(PF_reaction_writer(decomp_network_simple).write_PF_reactions() )

print('\n\n#### CORPSE-style network\n\n')
print(PF_reaction_writer(decomp_network_CORPSE).write_PF_reactions() )

import matplotlib.pyplot as pyplot
def make_nodecolors(nodes,POMcol='C0',microbecol='C1',DOMcol='C2',MAOMcol='C3'):
    colors=[]
    for node in nodes:
        if 'MICROBES' in node:
            colors.append(microbecol)
        elif 'MAOM' in node:
            colors.append(MAOMcol)
        elif 'DOM' in node:
            colors.append(DOMcol)
        else:
            colors.append(POMcol)
    return colors
    

pyplot.figure('Network comparison');pyplot.clf()
ax1=pyplot.subplot(121)
pos_simple=nx.layout.spectral_layout(decomp_network_simple)
nx.draw_networkx(decomp_network_simple,pos=pos_simple,with_labels=True,ax=ax1,nodes=decomp_network_simple.nodes,node_color=make_nodecolors(decomp_network_simple.nodes),arrowsize=20)
pyplot.title('Simple network')

ax2=pyplot.subplot(122)
pos_CORPSE=nx.layout.spectral_layout(decomp_network_CORPSE)
pos_CORPSE['MAOM'][0] += 0.5
pos_CORPSE['MAOM'][1] = pos_CORPSE['NECROMASS'][1]-0.3
nx.draw_networkx(decomp_network_CORPSE,pos=pos_CORPSE,with_labels=True,ax=ax2,nodes=decomp_network_CORPSE.nodes,node_color=make_nodecolors(decomp_network_CORPSE.nodes),arrowsize=20)
pyplot.title('CORPSE-style network')

pyplot.tight_layout()
pyplot.show()