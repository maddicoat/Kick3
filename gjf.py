from ase.atoms import Atoms
from ase.parallel import paropen


def write_gjf(fileobj, images):
    if isinstance(fileobj, str):
        fileobj = paropen(fileobj, 'w')

    if not isinstance(images, (list, tuple)):
        images = [images]

    symbols = images[0].get_chemical_symbols()
    natoms = len(symbols)
    for atoms in images:
        this_info = atoms.info
        fileobj.write('%5s%-20s \n' % ('%mem=',this_info["mem"]))
        fileobj.write('%7s \n' % ('%NoSave'))
        fileobj.write('%5s%8s%4s\n' % ('%chk=',this_info["filebase"],'.chk'))
        fileobj.write('%6s%-20s\n' % ('%nproc=',this_info["ncpus"]))
        fileobj.write('%13s \n\n' % ('#Method/basis'))
        fileobj.write('%8s \n\n' % (this_info["filebase"]))
        fileobj.write('%d %d\n' % (this_info["charge"], this_info["multiplicity"]))
        for s, (x, y, z) in zip(symbols, atoms.get_positions()):
            fileobj.write('%-2s %22.15f %22.15f %22.15f\n' % (s, x, y, z))
        if images[0].pbc.any():
            fileobj.write('%-5s %8.3f %8.3f %8.3f \n' %
                      ('Tv', images[0].get_cell()[0][0],
                       images[0].get_cell()[0][1],
                       images[0].get_cell()[0][2]))
            fileobj.write('%-5s %8.3f %8.3f %8.3f \n' %
                      ('Tv', images[0].get_cell()[1][0],
                       images[0].get_cell()[1][1],
                       images[0].get_cell()[1][2]))
            fileobj.write('%-5s %8.3f %8.3f %8.3f \n' %
                      ('Tv', images[0].get_cell()[2][0],
                       images[0].get_cell()[2][1],
                       images[0].get_cell()[2][2]))

