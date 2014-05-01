#!/usr/bin/python
from ase.atom import Atom
from ase.atoms import Atoms
from ase.calculators.neighborlist import *
from ase.data import *

from collections import namedtuple
from operator import attrgetter
from math import pi
from ase.io import read, write
import numpy as np
import argparse
import random
import operator

from conformer import random_conformer as rand_conf

def read_input(filename):
    """Method to read the Kick.in file and return the results in a dictionary"""
    
    really_big = 10000
    my_k3_params = {} 

    if isinstance(filename, str):
        myfile = open(filename)
    
    lines = myfile.readlines()     
       
    for line in lines:
        if (line.strip().startswith('#')):
            pass
        elif("=" in line):
            key, value = line.split("=")
            this_key = key.lower().strip()
            if("geom" in this_key):
                my_k3_params["num_geoms"] = int(value)
            elif("ncpu" in this_key or "ncore" in this_key):
                my_k3_params["ncpus"] = int(value)
            elif("mem" in this_key): #put this as a separate case in case we want to parse it later
                my_k3_params["mem"] = value.strip()
            elif("wall" in this_key or "time" in this_key):
                my_k3_params["wall"] = int(value)
            elif("cha" in this_key):
                my_k3_params["charge"] = int(value)
            elif("mult" in this_key):
                my_k3_params["multiplicity"] = int(value)
            elif("radius" in this_key):
                my_k3_params["radius"] = float(value)
            elif("boundary" in this_key):
                boundaries = [float(tmp) for tmp in value.split(',')]
                if(len(boundaries > 2)):
                    raise ValueError('More than two boundaries requested!')
                my_k3_params["boundaries"] = boundaries 
            elif("pbc" in this_key):
                cell = [float(tmp) for tmp in value.split(',')]
                pbc = [True for tmp in value.split(',')]
                if len(cell) > 3:
                    raise ValueError('Cell may be periodic in x,y,z dimensions and no more!')
                while len(cell) < 3:
                    cell.append(really_big)
                    pbc.append(False)
                my_k3_params["pbc"] = pbc
                my_k3_params["cell"] = cell
            elif("order" in this_key):
                if("AA" in value.strip()):
                    my_k3_params["order"] = False
                else:
                    my_k3_params["order"] = True
            elif("prog" in this_key):
                prog = value.lower().strip()
                if("gau" in prog or "g0" in prog):
                    my_k3_params["program"] = "gaussian"
                elif("dftb" in prog):
                    my_k3_params["program"] = "dftb"
                elif("xyz" in prog):
                    my_k3_params["program"] = "xyz"
                else:
                    raise ValueError('Only programs supported so far are gaussian and dftb+. \n Use "xyz" for output of simple coordinates')
            elif("save" in this_key):
                my_k3_params["saveall"] = True
                
    myfile.close()           
    #Now check for all essential stuff
    if("num_geoms" not in my_k3_params):
        raise ValueError('Number of geometries not set: NumGeom = ??')        
    #if("wall" not in my_k3_params):
    #    raise ValueError('per job walltime not defined: Time = ?? (in hours)') 
    if("radius" not in my_k3_params):
        raise ValueError('radius not defined: radius = ?? (in Angstrom)') 
    if("program" not in my_k3_params):
        raise ValueError('program not defined: supply either gaussian, dftb or xyz')
    if(my_k3_params["program"] == "gaussian"):
        if("charge" not in my_k3_params):
            raise ValueError('Charge not defined: Charge = ??')
        if("multiplicity" not in my_k3_params):
            raise ValueError('Multiplicity not defined: Mult = ??')
        if("ncpus" not in my_k3_params):
            raise ValueError('Number of cpus not set: NCPUs = ??')             
        if("mem" not in my_k3_params):
            raise ValueError('Memory not defined (supply units): Mem = ??')
    #my_k3_params[key.strip()] = value.strip()
            
   
    return my_k3_params

def read_system(filename):
    """Method to read the System.in file and set up the corresponding dictionaries"""

    SystemDict = namedtuple('SystemDict', 'number, frag_type')
    Tagpair = namedtuple('Tagpair', ['a', 'b']) #for flexible fragments
    system_def = {}
    system_master = {}
    system_tags = {}

    if isinstance(filename, str):
        myfile = open(filename)
    
    lines = myfile.readlines()

    #We have both atoms and fragments
    start_reading_atoms = False
    stop_reading_atoms = False
    start_reading_frags = False
    stop_reading_frags = False

    for line in lines:
        if (line.strip().startswith('#')):
            pass
        else:
            if ('Atoms' in line):
                start_reading_atoms = True
                continue
            if start_reading_atoms:
                if ('End' in line):
                    stop_reading_atoms = True
                    continue
            if(start_reading_atoms and not stop_reading_atoms):
                    atom, number_of = line.split()[:2]
                    system_def[atom] = SystemDict(int(number_of), "A")
            if('Fragments' in line):
                start_reading_frags = True
                continue
            if start_reading_frags:
                if ('End' in line):
                    stop_reading_atoms = True
                    continue
            if(start_reading_frags and not stop_reading_frags):
                    #print line
                    frag, number_of, frag_type = line.split()[:3]
                    print frag_type.upper()[0]
                    system_def[frag] = SystemDict(int(number_of), frag_type.upper()[0])
    myfile.close()

    #Now set up a dictionary with the actual atoms objects
    for species,info in system_def.items():
        if info.frag_type == "A":
            system_master[species] = Atoms('species', positions=[(0, 0, 0)])
        if info.frag_type == "R":
            filename = species + '.xyz'
            system_master[species] = read(filename)
        if info.frag_type == "L":
            filename = species + '.xyz'
            system_master[species] = read(filename)
        if info.frag_type == "F":
            tag_pairs = []
            filename = species + '.xyz'
            system_master[species] = read(filename)
            data_file = species + '.txt'
            print data_file
            try:
                with open(data_file) as f:
                    for line in f.readlines():
                        (first,second) =line.split()
                        this_pair = Tagpair(int(first.strip()),int(second.strip()))
                        tag_pairs.append(this_pair)
                    system_tags[species] = tag_pairs
            except IOError:
               raise IOError('File '+`data_file`+' for flexible fragment '+species+" doesn't exist")
    return system_def, system_master, system_tags

def check_coords(mol0,mol1):
    """Check distances between old and new molecules (Atoms objects) 
    We'll add the molecules together such that we can use the mic facility.
    mol0 is the 'current good' mol and carries the pbc information, if any."""
    
    keep_coords = True
    
    tmp_mol = mol0 + mol1
    mol0_natoms = len(mol0)
    #print "Check", len(mol0), len(mol1), len(tmp_mol), mol0_natoms
    for a0 in xrange(mol0_natoms-1):
        for a1 in xrange(mol0_natoms,len(tmp_mol)):
            #print a0,a1
            if tmp_mol[a0].symbol != "X" and tmp_mol[a1].symbol != "X":
                this_dist = tmp_mol.get_distance(a0, a1, mic = True)
                cutoff_dist = covalent_radii[tmp_mol[a0].number] + covalent_radii[tmp_mol[a1].number] + 0.3
                #print tmp_mol[a0].symbol, tmp_mol[a1].symbol,this_dist, cutoff_dist
                if this_dist < cutoff_dist:
                    keep_coords = False
                    #print "returning false, overlap"
                    return keep_coords

    #Second way to fail out if we have PBC (first way - i.e. overlap, is covered above by using mic=True)
    # The centroid of the addition must be within the box
    #The +/- 1/2 test used here assume we've already centred the cell
    com = mol1.get_center_of_mass()
    #print com
    #axes = [i for i in range(len(mol0.get_pbc())) if mol0.pbc[i]]
    for axis in xrange(3):
        if mol0.pbc[axis]:
            if com[axis] > mol0.cell[axis,axis] or com[axis] < -mol0.cell[axis,axis]:  #orthorhombic cell 
                keep_coords = False
                #print "returning false, PBC"
                return keep_coords

    #print "keeping coords ", keep_coords            
    return keep_coords
        
    
    

#def get_intermolecular_distance(mol0, mol1, a0, a1):
#    """Return distance between two atoms in different molecules."""
#
#    R0 = mol0.arrays['positions']
#    R1 = mol1.arrays['positions']
#    D = R1[a1] - R0[a0]
#    return np.linalg.norm(D)

#------------------------------------------------------------------------------

system_input = "System.in"
main_input = "Kick.in"
restart_file = "Restart.data"

try:
    with open(restart_file) as f:
        line = f.readlines(1)
        start = int(line[0].strip()) + 1
except IOError:
   start = 0


system_def, system_master, system_tags = read_system(system_input)
k3_params = read_input(main_input)

print k3_params
print system_def
print system_master
print system_tags

#Get total number of things
total_objects = sum([v.number for v in system_def.values()])
print "TO= ", total_objects

#create a list of fragments in the order we'll add them
frag_list = []
tmp_sysdef ={} #We can loop destructively over this
for k,v in system_def.items():
    tmp_sysdef[k] = int(v.number)
    

print tmp_sysdef
if k3_params["order"]:
#while any (numleft > 0.0 for numleft in tmp_sysdef.values()): #works on command line, not in Spyder
    while sum(tmp_sysdef.values()) != 0:
        #print "In while loop", tmp_sysdef.values()
        for k,v in sorted(tmp_sysdef.iteritems(), key=operator.itemgetter(1)):
            if v > 0:
                #print "In inner if"
                frag_list.append(k)
                tmp_sysdef[k] -= 1
                #print tmp_sysdef.values()
else:
    for k,v in sorted(tmp_sysdef.iteritems(), key=operator.itemgetter(1)):
        for tmp in range(v):
            frag_list.append(k)            
    
print frag_list    
######
for ci in xrange(k3_params["num_geoms"]):
    current = start + ci #LOOP over numgeom GOES HERE
    this_file = "%04d" % current
    if k3_params["program"] == 'dftb':
        filename='Kick'+this_file+'.gen'  
    elif k3_params["program"] == 'gaussian':
        filename='Kick'+this_file+'.gjf'
        filebase='Kick'+this_file
    else:
        filename='Kick'+this_file+'.xyz'  

    #Now, create a new Atoms object
    #place a dummy at the origin and centre it (takes care of pbc)
    if "pbc" in k3_params.keys():
        new_mol = Atoms(pbc=k3_params["pbc"], cell=k3_params["cell"])
        dummy_mol = Atoms(pbc=k3_params["pbc"], cell=k3_params["cell"])
    else:
        new_mol = Atoms()
        dummy_mol = Atoms()
        
    #special processing for gaussian
    if k3_params["program"] == 'gaussian':
        new_mol.info["mem"] = k3_params["mem"]
        new_mol.info["ncpus"] = k3_params["ncpus"]
        new_mol.info["filebase"] = filebase
        new_mol.info["charge"] = k3_params["charge"]
        new_mol.info["multiplicity"] = k3_params["multiplicity"]
    
    origin = Atom('X', position=(0,0,0))
    
    new_mol += origin
    dummy_mol += origin
    
    axes = [i for i in range(len(new_mol.get_pbc())) if new_mol.pbc[i] ]
    
    new_mol.center(axis=axes)
    dummy_mol.center(axis=axes)
    new_mol.pop()
    
    #Now loop over the list
    for this_frag in frag_list:
        #print "Now here!"
        keep_coords = False
        attempts = 0
        if system_def[this_frag].frag_type == "L": #Locked fragment, dump it no questions asked
            new_mol += system_master[this_frag]
            this_com = system_master[this_frag].get_center_of_mass()
            dummy_mol += this_com
        else:
            #choose a centroid
            #print dummy_mol
            this_origin = random.choice(dummy_mol).position
            #print this_origin
            attempts = 0
            while not keep_coords:            
                #translation = this_origin + (np.random.random_sample(3) * k3_params["radius"] * 2 - k3_params["radius"])
                translation = np.random.random_sample(3) * k3_params["radius"] * 2 - k3_params["radius"]          
                translation += this_origin
                #print "new_origin = ",translation
                this_phi = random.uniform(-pi,pi)
                this_theta = random.uniform(-pi,pi)
                this_psi = random.uniform(-pi,pi)
                tmp_mol = system_master[this_frag].copy()
                if system_def[this_frag].frag_type == "F":
                    #print this_frag, " is flexible"
                    rand_conf(tmp_mol,system_tags[this_frag])
                if system_def[this_frag].frag_type != "A":
                    tmp_mol.rotate_euler(center='COM',phi=this_phi,theta=this_theta,psi=this_psi)
                tmp_com = tmp_mol.get_center_of_mass()
                tmp_mol.translate(translation - tmp_com)    
                attempts+=1
                keep_coords = check_coords(new_mol,tmp_mol)
                #print keep_coords
                if keep_coords:
                    #print "Here!"
                    new_mol+=tmp_mol
                    tmp_com = tmp_mol.get_center_of_mass()
                    dummy_mol += Atom('X',position=tmp_com)
                    break
                elif attempts == 5:
                    #Select a new origin
                    #This is in case we've selected something near the centre, leaving no legal placement 
                    #for the new fragment
                    this_origin = random.choice(dummy_mol).position
                    attempts = 0
                
    write(filename,new_mol)
    #write('centroids.xyz',dummy_mol)

#dump for restart
with open(restart_file, 'w') as f:
    f.write('%d\n' % (current))