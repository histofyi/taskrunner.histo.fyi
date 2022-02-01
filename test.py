from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout


pdb = PDBFile('input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
print('Forcefield applied')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
print('System created')
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
print('Simulation ready to start')
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
print('Minimization complete')
simulation.reporters.append(PDBReporter('output.pdb', 100))
simulation.reporters.append(StateDataReporter(stdout, 100, step=True, potentialEnergy=True, temperature=True))
simulation.step(10000)
