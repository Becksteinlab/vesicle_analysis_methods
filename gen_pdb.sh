module load gromacs/4.6.6
cd systems
for d in */ ; do
    if [ -e $d/emin/emin.pdb ]
        then
            rm $d/emin/emin.pdb
    fi

    editconf -f $d/emin/emin.gro -o $d/emin/emin.pdb
done
