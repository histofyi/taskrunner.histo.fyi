from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

import datetime
import json

statistics = {
    'items': []
}


steps = 100000
report_frequency = 1000
inputfile = 'input.pdb'
outputfile = 'output/test.pdb'


"""
This function adds a row to the item array in the statistics dictionary

Args:
    item (string) : the name of the activity which has just been completed
"""
def add_stats_row(item: str):
    statistics['items'].append({'step':item, 'started_at':datetime.datetime.now()})
    if item == 'start' or item == 'end':
        statistics[item] = datetime.datetime.now()


"""
This function processes the statistics dictionary at the end of a simulation

Args:
    steps (integer) : the number of steps in the simulation
    inputfile (string) : the filename of the structure being simulated
    outputfile (string) : the filename of the structure after simulation
"""
def process_statistics(steps : int, inputfile : str, outputfile : str):
    statistics['steps'] = steps
    statistics['input'] = inputfile
    statistics['output'] = outputfile
    statistics['duration'] = (statistics['end'] - statistics['start']).total_seconds()
    statistics['start'] = statistics['start'].isoformat()
    statistics['end'] = statistics['end'].isoformat()
    for item in statistics['items']:
        item['started_at'] = item['started_at'].isoformat()
    f = open('output/logfile.json', 'w')
    f.write(json.dumps(statistics, indent=4, sort_keys=True))
    f.close()

    print(statistics)


"""
This function runs a simulation using OpenMM

Args:
    steps (integer) : the number of steps in the simulation
    inputfile (string) : the filename of the structure being simulated
    outputfile (string) : the filename of the structure after simulation
"""
# adapted from the OpenMM simulatePDB example
def run_simulation(steps : int, inputfile : str, outputfile : str, report_frequency : int =1000):
    add_stats_row('start')
    pdb = PDBFile('input.pdb')
    add_stats_row('force_field')
    print('Forcefield applied')
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
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
    simulation.reporters.append(PDBReporter(outputfile, report_frequency))
    simulation.reporters.append(StateDataReporter(stdout, report_frequency, step=True, potentialEnergy=True, temperature=True))
    simulation.step(steps)
    add_stats_row('end')
    process_statistics(steps, inputfile, outputfile)


# run the simulation
run_simulation(steps, inputfile, outputfile, report_frequency=report_frequency)