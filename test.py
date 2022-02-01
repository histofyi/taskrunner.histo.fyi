from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout

import datetime

statistics = {}


def add_stats_row(item):
    statistics[item] = datetime.datetime.now().isoformat()



add_stats_row('start')
pdb = PDBFile('input.pdb')
add_stats_row('pdb_loaded')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
add_stats_row('force_field')
print('Forcefield applied')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
add_stats_row('system')
print('System created')
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
add_stats_row('simulation_ready')
print('Simulation ready to start')
simulation.context.setPositions(pdb.positions)
add_stats_row('set_positions')
simulation.minimizeEnergy()
add_stats_row('minimisation_complete')
print('Minimization complete')
simulation.reporters.append(PDBReporter('output.pdb', 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
add_stats_row('end')

print(statistics)
