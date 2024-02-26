import sys 
import numpy as np 
import subprocess,os  
print("You are running the program:", sys.argv[0], " to geneate lammps input files" )
#ntimes = [2,2] 
boxa = 40.0 
d    = 1.5  # using in pakckmol to get box size reduction  
boxb = boxa 
boxc = boxa 
ntimes = sys.argv[1::2]   # even elements are numbers of molecules
ntimes = [int(i) for i in ntimes] 
files  = sys.argv[::2]  # odd elements are molecule files
del files[0]  
nfiles  = len(files) 
print(ntimes)
print(files) 
print(nfiles) 
# number of fragments of each molecule 2 benz, 2pc , 3 act etc.. 
# order of the files should match 
print("Number of lammps files to read ", nfiles) 
#print("Argument List:", str(sys.argv[1:])) 
f = open("final.lmp", "w") 
f.write('# created by simil\n\n')
f.write('log final.log append\n')
f.write('units real\n')
f.write('atom_style full\n')
f.write('boundary p p p\n\n')

f.write('pair_style lj/cut 3.0\n')
f.write('kspace_style none\n')
f.write('dielectric 1.000000\n')
f.write('bond_style harmonic\n')
f.write('angle_style harmonic\n')
f.write('dihedral_style fourier\n')   # change dihedral_style  
f.write('improper_style cvff\n')
f.write('special_bonds amber\n\n')
f.write('pair_modify mix arithmetic\n')
f.write('neighbor 2.0 bin\n')
f.write('neigh_modify delay 0 every 1 check yes\n')
f.write('read_data final_data.lmp\n')
f.write('# read_restart restart1.lmp\n\n')

flist = [] 
occ = 0 
for i in files:  
#   print(i) 
    f1 = open(i, "r")
    txt = f1.read() 
    x = txt.count("pair_coeff") 
    f2 = txt.splitlines()  
    f1.close() 
    for line in f2: 
        if "read_data" in line: 
           f1 = line.split() 
           flist.append(f1[1]) 
        if "Pair_coeff" in line: 
           occ = occ + 1 
           y = line.split() 
           y[1] = str(occ)  
           y[2] = str(occ)  
           jo = " ".join(y)
#          print(jo) 
           f.write(jo) 
           f.write('\n') 
#print(nocc) 
f.write('min_style cg\n')
f.write('minimize 0.000100 0.000001 10000 100000  \n')
f.write('reset_timestep 0 \n')
f.write('dump dump0 all custom 1000 eq1.dump id type mol x y z ix iy iz vx vy vz  \n')
f.write('dump xtc0 all xtc 1000 eq1.xtc  \n')
f.write('\n')

f.write('# nvt  \n')
f.write('timestep 1.0\n\n')
f.write('fix md1 all nvt temp 300.000000 300.000000 100.000000  \n')
f.write('run 1000000  \n')
f.write('unfix md1 \n')

f.write('thermo_style custom step time temp press enthalpy etotal ke pe ebond eangle edihed eimp evdwl ecoul elong etail vol lx ly lz density pxx pyy pzz pxy pxz pyz \n')
f.write('thermo_modify flush yes \n')
f.write('dump TRAJ all custom 1000 dump.lammpstrj id mol type element q xu yu zu\n')
#f.write('dump_modify TRAJ element')
f.write('write_data data.eq.lmp\n')
f.write('write_dump all custom eq1_last.dump id x y z xu yu zu vx vy vz fx fy fz modify sort id \n')
f.write('write_data eq2.data\n')

f.close() 
print(flist) 
oc = 0 
atoms       = [] 
bonds       = [] 
angles      = [] 
dihedrals   = [] 
improper    = [] 
atom_types  = [] 
bond_types  = [] 
angle_types = [] 
dihedral_types = [] 
improper_types = [] 
for i in flist: 
    f1 = open(i, "r") 
#   print(i) 
    impfind = 0 
    for line in f1: 
        if "atoms" in line: 
            rline =  line.split() 
            atoms.append(int(rline[0]))
        if "bonds" in line: 
            rline =  line.split() 
            bonds.append(int(rline[0]))
        if "angles" in line: 
            rline =  line.split() 
            angles.append(int(rline[0]))
        if "dihedrals" in line: 
            rline =  line.split() 
            dihedrals.append(int(rline[0]))
        if "impropers" in line: 
            rline =  line.split() 
            improper.append(int(rline[0]))
        if "atom types" in line: 
            rline =  line.split() 
            atom_types.append(int(rline[0]))
        if "bond types" in line: 
            rline =  line.split() 
            bond_types.append(int(rline[0]))
        if "angle types" in line: 
            rline =  line.split() 
            angle_types.append(int(rline[0]))
        if "dihedral types" in line: 
            rline =  line.split() 
            dihedral_types.append(int(rline[0])) 
        if "improper types" in line: 
            impfind = 1 
            rline =  line.split() 
            improper_types.append(int(rline[0])) 
    if impfind==0: 
         improper_types.append(0) 
         improper.append(0)
    print("imp",improper)
#print(atoms,ntimes)  
#f.write('\n')
satoms = sum(np.multiply(atoms,ntimes))  
sbonds = sum(np.multiply(bonds,ntimes))
sangles = sum(np.multiply(angles,ntimes)) 
sdihedrals = sum(np.multiply(dihedrals,ntimes)) 
simpropar = sum(np.multiply(improper,ntimes)) 
print(dihedrals,ntimes) 
print('improper',improper,ntimes) 
satom_types = sum(atom_types) 
sbond_types = sum(bond_types) 
sangle_types = sum(angle_types) 
sdihedral_types = sum(dihedral_types) 
simproper_types = sum(improper_types) 
print('dihedral_types',dihedral_types)
print('improper_types',improper_types)
f = open('final_data.lmp', "w") 
f.write('created by simil\n\n')
f.write('{:d} atoms\n'.format(satoms)) 
f.write('{:d} bonds\n'.format(sbonds))
f.write('{:d} angles\n'.format(sangles))
f.write('{:d} dihedrals\n'.format(sdihedrals))
f.write('{:d} impropers\n'.format(simpropar))
f.write('{:d} atom types\n'.format(satom_types))
f.write('{:d} bond types\n'.format(sbond_types))
f.write('{:d} angle types\n'.format(sangle_types))
f.write('{:d} dihedral types\n'.format(sdihedral_types)) 
f.write('{:d} improper types\n'.format(simproper_types)) 

f.write('{:f} {:f} xlo xhi\n'.format(0.0,boxa)) 
f.write('{:f} {:f} ylo yhi\n'.format(0.0,boxb))
f.write('{:f} {:f} zlo zhi\n'.format(0.0,boxc)) 

print(satoms, "atoms") 
print(sbonds,"bonds" ) 
print(sangles,"angles") 
print(sdihedrals,"dihedrals" ) 
print(satom_types,"atom types" ) 
print(sbond_types ,"bond types") 
print(sangle_types ,"angle types") 
print(sdihedral_types,"dihedral types" ) 
print(simproper_types,"improper types" ) 
# reading masses over 
nfile = 0 
imass = 1 
gentype  = [] 
#print("\n")
f.write('\nMasses\n\n')
#print("Masses \n")
for i in flist: 
    gentype1 = [] 
    f1 = open(i, "r") 
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0  
    for k in range(nlines): 
        line = f1.readline()
        if "Masses" in line: 
           line = f1.readline()
           ifound = 1
           for i in range(atom_types[nfile]): 
               line = f1.readline().split() 
               line0   = line.copy() 
               gentype1.append(line0)
               line[0] = str(imass) 
               jo = " ".join(line)
#              print(jo) 
               f.write(jo) 
               f.write("\n") 
               imass = imass + 1 
        elif ifound==1: 
             break 
    nfile = nfile + 1       
    f1.close() 
    gentype.append(gentype1) 
#print(gentype) 
# bond 
# Pair Coeffs
f.write('\nPair Coeffs\n\n')
nfile = 0 
ivdw = 1 
gentype  = [] 
for i in flist:
    gentype1 = []
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Pair Coeffs" in line:
           line = f1.readline()
           ifound = 1
           for i in range(atom_types[nfile]):
               line = f1.readline().split()
               line0   = line.copy()
               gentype1.append(line0)
               line[0] = str(ivdw)
               jo = " ".join(line)
               f.write(jo)
               f.write("\n")
               ivdw = ivdw + 1
        elif ifound==1:
             break
    nfile = nfile + 1
    f1.close()
    gentype.append(gentype1)
#print("\n")
f.write('\nBond Coeffs\n\n')
#print("Bond Coeffs \n")
nfile = 0
imass = 1
for i in flist:
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Bond Coeffs" in line:
           line = f1.readline()
           ifound = 1
           for i in range(bond_types[nfile]):
               line = f1.readline().split()
               line[0] = str(imass)
               jo = " ".join(line)
#              print(jo)
               f.write(jo)
               f.write("\n")
               imass = imass + 1
        elif ifound==1:
             break
    nfile = nfile + 1
    f1.close() 
#Angle coefficient 
#print("\n")
f.write('\nAngle Coeffs\n\n')
#print("Angle \n")
nfile = 0
imass = 1
for i in flist:
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Angle Coeffs" in line:
           line = f1.readline()
           ifound = 1
           for i in range(angle_types[nfile]):
               line = f1.readline().split()
               line[0] = str(imass)
               jo = " ".join(line)
#              print(jo)
               f.write(jo)
               f.write("\n")
               imass = imass + 1
        elif ifound==1:
             break
    nfile = nfile + 1
    f1.close() 
# Dihedral Coeffs 
#print("\n")
f.write('\nDihedral Coeffs\n\n')
#print("Dihedral  \n")
nfile = 0
imass = 1
for i in flist:
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Dihedral Coeffs" in line:
           line = f1.readline()
           ifound = 1
           for i in range(dihedral_types[nfile]):
               line = f1.readline().split()
               line[0] = str(imass)
               jo = " ".join(line)
               f.write(jo)
               f.write("\n")
               imass = imass + 1
        elif ifound==1:
             break
    nfile = nfile + 1
    f1.close() 
# Improper Coeffs  
f.write('\nImproper Coeffs\n\n')
nfile = 0
imass = 1
for i in flist:
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Improper Coeffs" in line:
           line = f1.readline()
           ifound = 1
           for i in range(improper_types[nfile]):
               line = f1.readline().split()
               line[0] = str(imass)
               jo = " ".join(line)
               f.write(jo)
               f.write("\n")
               imass = imass + 1
        elif ifound==1:
             break
    nfile = nfile + 1
    f1.close() 
# NOw runpack mol 
nfile = 0
fpack = [] 
nmol  = [] 
for i in flist:
#   print(i)
    fn = i.split(".")[0]+".xyz" 
    fpack.append(fn) 
    print(fn) 
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    a = [] 
    x = [] 
    y = [] 
    z = [] 
#   print(gentype[nfile]) 
    for k in range(nlines):
        line = f1.readline()
        if "Atoms" in line:
           ifound = 1
           line = f1.readline()
#          print(atoms[nfile])
           for i in range(atoms[nfile]):
               line = f1.readline().split()
#              print(line) 
               tat  = line[2]  
               prn = 0 
               for kk in gentype[nfile]: 
#                  print(kk) 
                   if kk[0]==tat:  
                       atomget = kk[-1].split(",")[0]  
                       if atomget == "na": 
                           atomget= "Na"
                       else: 
                           atomget =  atomget[0].upper() 
                   #   print(atomget) 
                       prn = atomget # kk[3][0]  
#              print(prn,float(line[4]),float(line[5]),float(line[6])) 
               a.append(prn) 
               x.append(float(line[4])) 
               y.append(float(line[5])) 
               z.append(float(line[6])) 
        elif ifound==1:
             break
    nfile = nfile + 1
    f1.close() 
    fc = open(fn, "w")
    fc.write(str(len(a)))
    fc.write("\n")
    fc.write("\n")
    nmol.append(len(a)) 
    for i in range(len(a)):
        fc.write('{:5s} {:15.6f} {:15.6f} {:15.6f}\n'.format(a[i], x[i], y[i], z[i]))
#       print('{:5s} {:15.6f} {:15.6f} {:15.6f}\n'.format(a[i], x[i], y[i], z[i]))
    fc.close() 

print(fpack) 
print(nmol) 
fp= open("packmol.inp", "w") 
fp.write('# created by simil\n')
tol=2.5 
fp.write(f'tolerance {tol:3.1f}\n')
fp.write('filetype xyz\n')
fp.write(f'output simbox.xyz \n')
#fp.write('\n')

for i in range(len(fpack)): 
    nm1 = ntimes[i]  
    xyzfile = fpack[i] 
    fp.write(f'\nstructure {xyzfile}\n') 
    fp.write(f'  number {nm1}\n') 
    fp.write('  inside box {:.4f} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n'.format(d, d, d, boxa - d, boxb - d, boxc - d))
    fp.write('end structure\n')
fp.close() 
# end packmol generation 
# packmol input generated : now run the packmol command 
print("running packmol") 
os.system("cat packmol.inp") 
#completed_process = subprocess.Popen(['/home/thoms0a/Downloads/packmol/packmol-20.14.0/packmol  < packmol.inp']) #capture_output=True, text=True, shell = True)
os.system("packmol  < packmol.inp > pack.out") 
#os.system("/home/thoms0a/Downloads/packmol/packmol-20.14.0/packmol  < packmol.inp") 
# opening simbox.xyz 
fm= open("simbox.xyz", "r") 
line = fm.readline().split()  
natoms = int(line[0])
print('from simbox.out #of atoms',natoms) 
line = fm.readline()
atomtype = [] 
x = [] 
y = [] 
z = [] 
for i in range(natoms): 
    line = fm.readline().split()  
#   print(line,i) 
#   atomtype[i]= line[0] 
    x.append(line[1]) 
    y.append(line[2]) 
    z.append(line[3]) 
fm.close() 
# Atoms coordinates 
#print("\n")
#print("Atoms \n")
f.write('\nAtoms\n\n')
nfile = 0
imass = 1
iadd  = 1 
iatom_types = 0 
for i in flist:
    f1 = open(i, "r")
    f1.seek(0) 
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
#   print(nlines,i) 
    for k in range(nlines):
        line = f1.readline()
#       print(line) 
        if "Atoms" in line:
           line = f1.readline()
#          print(line) 
           ifound = 1
           life = [] 
           for i in range(atoms[nfile]):
#              print(line) 
               line = f1.readline().split()
               life.append(line)
        # add atom type seperately : else it will add many times 
           for line in life:  
               line[2] = str((int(line[2]) + iatom_types)) 
           for m in range(ntimes[nfile]): 
               for line in life:  
                   line[0] = str(imass)
                   line[1] = str(iadd)
#                  print(imass) 
#                  line[4] = x[imass]
#                  line[5] = y[imass]
#                  line[6] = z[imass]
#                  print(line) 
                   jo = " ".join(line)
#                  print(jo)
                   f.write(jo)
                   f.write("\n")
                   imass = imass + 1
               iadd = iadd + 1
        elif ifound==1:
             break
    iatom_types = iatom_types + atom_types[nfile] 
    nfile = nfile + 1
    f1.close() 
#print("\n")

#print("Bonds \n")
f.write('\nBonds\n\n')
ibond_types = 0 
ibond = 1
nfile = 0
ibonda = 0 
for i in flist:
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
#   print(nlines) 
    for k in range(nlines):
        line = f1.readline()
        if "Bonds" in line:
           line = f1.readline()
           ifound = 1
           life = []
           for i in range(bonds[nfile]):
               line = f1.readline().split()
               life.append(line)
        # add bomd type seperately : else it will add many times
        #  print(life)  
           for line in life:
               line[1] = str((int(line[1]) + ibond_types))
           for m in range(ntimes[nfile]):
               for line in life:
                   line0   = line.copy() 
                   line0[0] = str(ibond)
        #          print(ibonda,line[2],line[3]) 
                   line0[2] = str( int(line[2]) + ibonda )
                   line0[3] = str( int(line[3]) + ibonda )
                   jo = " ".join(line0)
#                  print(jo)
                   f.write(jo)
                   f.write("\n")
                   ibond = ibond + 1
               ibonda  = ibonda +  atoms[nfile] 
#              iadd = iadd + 1
        elif ifound==1:
             break
    ibond_types = ibond_types + bond_types[nfile]
    nfile = nfile + 1
    f1.close() 

#print("\n")
#print("Angles \n")
f.write('\nAngles\n\n') 
iangle_types = 0
iangle= 1
nfile = 0
ianglea = 0
for i in flist:
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Angles" in line:
           line = f1.readline()
           ifound = 1
           life = []
           for i in range(angles[nfile]):
               line = f1.readline().split()
               life.append(line)
        # add bomd type seperately : else it will add many times
        #  print(life)
           for line in life:
               line[1] = str((int(line[1]) + iangle_types))
           for m in range(ntimes[nfile]):
               for line in life:
                   line0   = line.copy()
                   line0[0] = str(iangle)
        #          print(ibonda,line[2],line[3])
                   line0[2] = str( int(line[2]) + ianglea )
                   line0[3] = str( int(line[3]) + ianglea )
                   line0[4] = str( int(line[4]) + ianglea )
                   jo = " ".join(line0)
#                  print(jo)
                   f.write(jo)
                   f.write("\n")
                   iangle = iangle + 1
               ianglea  = ianglea +  atoms[nfile]
#              iadd = iadd + 1
        elif ifound==1:
             break
    iangle_types = iangle_types + angle_types[nfile]
    nfile = nfile + 1
    f1.close() 

#print("\n")
#print("Dihedrals \n")
f.write('\nDihedrals\n\n')
idihedral_types = 0
idihed = 1
nfile = 0
idiheda = 0
for i in flist:
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Dihedrals" in line:
           line = f1.readline()
           ifound = 1
           life = []
           for i in range(dihedrals[nfile]):
               line = f1.readline().split()
               life.append(line)
        # add bomd type seperately : else it will add many times
        #  print(life)
           for line in life:
#              print(line[1],idihedral_types) 
               line[1] = str((int(line[1]) + idihedral_types))
#              print(line[1],idihedral_types) 
           for m in range(ntimes[nfile]):
               for line in life:
                   line0   = line.copy()
                   line0[0] = str(idihed)
                   line0[2] = str( int(line[2]) + idiheda )
                   line0[3] = str( int(line[3]) + idiheda )
                   line0[4] = str( int(line[4]) + idiheda )
                   line0[5] = str( int(line[5]) + idiheda )
                   jo = " ".join(line0)
#                  print(jo)
                   f.write(jo)
                   f.write("\n")
                   idihed = idihed + 1
               idiheda  = idiheda +  atoms[nfile]
#              iadd = iadd + 1
        elif ifound==1:
             break
    idihedral_types = idihedral_types + dihedral_types[nfile]
    nfile = nfile + 1
    f1.close() 

#print("\n")
#print("Dihedrals \n")
f.write('\nIMPROPER\n\n')
idihedral_types = 0
idihed = 1
nfile = 0
idiheda = 0
print("testing improper")
for i in flist:
#   print(i) 
    f1 = open(i, "r")
    with open(i, 'r') as fp:
       nlines = len(fp.readlines())
    ifound = 0
    for k in range(nlines):
        line = f1.readline()
        if "Impropers" in line:
           ifound = 1
           line = f1.readline()
    #      print(line) 
           life = []
           for i in range(improper[nfile]):
               line = f1.readline().split()
               life.append(line)
        # add bond type seperately : else it will add many times
    #      print('life',life)
           for line in life:
#              print(line[1],idihedral_types) 
               line[1] = str((int(line[1]) + idihedral_types))
#              print(line[1],idihedral_types) 
           for m in range(ntimes[nfile]):
               for line in life:
                   line0   = line.copy()
                   line0[0] = str(idihed)
                   line0[2] = str( int(line[2]) + idiheda )
                   line0[3] = str( int(line[3]) + idiheda )
                   line0[4] = str( int(line[4]) + idiheda )
                   line0[5] = str( int(line[5]) + idiheda )
                   jo = " ".join(line0)
                   print(jo)
                   f.write(jo)
                   f.write("\n")
                   idihed = idihed + 1
               idiheda  = idiheda +  atoms[nfile]
#              iadd = iadd + 1
        elif ifound==1:
             break
    idihedral_types = idihedral_types + dihedral_types[nfile]
    nfile = nfile + 1
    f1.close() 
