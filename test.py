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


def add_stats_row(item):
    statistics['items'].append({'step':item, 'started_at':datetime.datetime.now()})
    if item == 'start' or item == 'end':
        statistics[item] = datetime.datetime.now()


def process_statistics(steps, inputfile, outputfile):
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


def run_simulation(steps, inputfile, outputfile, report_frequency=1000):
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

run_simulation(steps, inputfile, outputfile, report_frequency=report_frequency)