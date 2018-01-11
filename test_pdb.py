from fragment_capping.helpers.molecule import molecule_from_pdb_str, Molecule, Atom

PDBS = {
    'ethanol': '''HEADER    UNCLASSIFIED                            21-Sep-17
TITLE     ALL ATOM STRUCTURE FOR MOLECULE LIG                                   
AUTHOR    GROMOS AUTOMATIC TOPOLOGY BUILDER REVISION 2017-09-18 11:12:19
AUTHOR   2  http://compbio.biosci.uq.edu.au/atb
HETATM    1   H1 X5VY    0      -1.341  -0.944   0.839  1.00  0.00           H
HETATM    2   C1 X5VY    0      -1.282  -0.256  -0.015  1.00  0.00           C
HETATM    3   H3 X5VY    0      -2.159   0.402   0.024  1.00  0.00           H
HETATM    4   H2 X5VY    0      -1.332  -0.847  -0.936  1.00  0.00           H
HETATM    5   C2 X5VY    0       0.004   0.563   0.032  1.00  0.00           C
HETATM    6   H5 X5VY    0       0.052   1.246  -0.824  1.00  0.00           H
HETATM    7   H4 X5VY    0       0.028   1.180   0.943  1.00  0.00           H
HETATM    8   O1 X5VY    0       1.182  -0.244  -0.061  1.00  0.00           O
HETATM    9   H6 X5VY    0       1.199  -0.816   0.724  1.00  0.00           H
CONECT    1    2
CONECT    2    1    3    4    5
CONECT    3    2
CONECT    4    2
CONECT    5    2    6    7    8
CONECT    6    5
CONECT    7    5
CONECT    8    5    9
CONECT    9    8
END''',
    'CHEMBL485374':  '''HEADER    UNCLASSIFIED                            21-Sep-17
TITLE     ALL ATOM STRUCTURE FOR MOLECULE UNL                                   
AUTHOR    GROMOS AUTOMATIC TOPOLOGY BUILDER REVISION 2017-09-18 11:12:19
AUTHOR   2  http://compbio.biosci.uq.edu.au/atb
HETATM    1   H1 DFWA    0       7.498  -0.961  -0.141  1.00  0.00           H
HETATM    2   C1 DFWA    0       6.425  -1.166  -0.104  1.00  0.00           C
HETATM    3   H3 DFWA    0       6.238  -1.814   0.760  1.00  0.00           H
HETATM    4   H2 DFWA    0       6.144  -1.721  -1.013  1.00  0.00           H
HETATM    5   N1 DFWA    0       5.719   0.097   0.038  1.00  0.00           N
HETATM    6   H4 DFWA    0       6.182   0.888  -0.390  1.00  0.00           H
HETATM    7   C2 DFWA    0       4.336   0.170   0.026  1.00  0.00           C
HETATM    8  C15 DFWA    0       3.517  -0.966   0.207  1.00  0.00           C
HETATM    9  H16 DFWA    0       3.965  -1.946   0.331  1.00  0.00           H
HETATM   10  C14 DFWA    0       2.131  -0.849   0.221  1.00  0.00           C
HETATM   11  H15 DFWA    0       1.542  -1.752   0.364  1.00  0.00           H
HETATM   12   C5 DFWA    0       1.483   0.392   0.058  1.00  0.00           C
HETATM   13   C6 DFWA    0       0.030   0.562   0.040  1.00  0.00           C
HETATM   14   H7 DFWA    0      -0.298   1.597  -0.064  1.00  0.00           H
HETATM   15   C7 DFWA    0      -0.911  -0.407   0.115  1.00  0.00           C
HETATM   16   H8 DFWA    0      -0.587  -1.446   0.178  1.00  0.00           H
HETATM   17   C8 DFWA    0      -2.365  -0.233   0.083  1.00  0.00           C
HETATM   18  C13 DFWA    0      -3.192  -1.361  -0.088  1.00  0.00           C
HETATM   19  H14 DFWA    0      -2.732  -2.343  -0.181  1.00  0.00           H
HETATM   20  C12 DFWA    0      -4.577  -1.258  -0.162  1.00  0.00           C
HETATM   21  H13 DFWA    0      -5.182  -2.152  -0.301  1.00  0.00           H
HETATM   22  C11 DFWA    0      -5.208  -0.006  -0.047  1.00  0.00           C
HETATM   23   N2 DFWA    0      -6.600   0.104  -0.043  1.00  0.00           N
HETATM   24  H12 DFWA    0      -7.086  -0.644  -0.526  1.00  0.00           H
HETATM   25  H11 DFWA    0      -6.954   1.009  -0.334  1.00  0.00           H
HETATM   26  C10 DFWA    0      -4.394   1.130   0.138  1.00  0.00           C
HETATM   27  H10 DFWA    0      -4.860   2.108   0.241  1.00  0.00           H
HETATM   28   C9 DFWA    0      -3.011   1.015   0.208  1.00  0.00           C
HETATM   29   H9 DFWA    0      -2.423   1.916   0.364  1.00  0.00           H
HETATM   30   C4 DFWA    0       2.315   1.519  -0.119  1.00  0.00           C
HETATM   31   H6 DFWA    0       1.856   2.498  -0.252  1.00  0.00           H
HETATM   32   C3 DFWA    0       3.698   1.420  -0.132  1.00  0.00           C
HETATM   33   H5 DFWA    0       4.305   2.312  -0.272  1.00  0.00           H
CONECT    1    2
CONECT    2    1    3    4    5
CONECT    3    2
CONECT    4    2
CONECT    5    2    6    7
CONECT    6    5
CONECT    7    5    8   32
CONECT    8    7    9   10
CONECT    9    8
CONECT   10    8   11   12
CONECT   11   10
CONECT   12   10   13   30
CONECT   13   12   14   15
CONECT   14   13
CONECT   15   13   16   17
CONECT   16   15
CONECT   17   15   18   28
CONECT   18   17   19   20
CONECT   19   18
CONECT   20   18   21   22
CONECT   21   20
CONECT   22   20   23   26
CONECT   23   22   24   25
CONECT   24   23
CONECT   25   23
CONECT   26   22   27   28
CONECT   27   26
CONECT   28   17   26   29
CONECT   29   28
CONECT   30   12   31   32
CONECT   31   30
CONECT   32    7   30   33
CONECT   33   32
END''',
    'warfarin': '''HEADER    UNCLASSIFIED                            31-Aug-17
TITLE     ALL ATOM STRUCTURE FOR MOLECULE WR0                                   
AUTHOR    GROMOS AUTOMATIC TOPOLOGY BUILDER REVISION 2017-07-03 14:53:07
AUTHOR   2  http://compbio.biosci.uq.edu.au/atb
HETATM    1  H16 AOOI    0       1.659  -3.906   2.643  1.00  0.00           H
HETATM    2  C14 AOOI    0       1.470  -3.838   1.570  1.00  0.00           C
HETATM    3  H14 AOOI    0       2.216  -4.415   1.013  1.00  0.00           H
HETATM    4  H15 AOOI    0       0.490  -4.283   1.351  1.00  0.00           H
HETATM    5  C13 AOOI    0       1.449  -2.394   1.132  1.00  0.00           C
HETATM    6   O4 AOOI    0       1.357  -1.493   1.963  1.00  0.00           O
HETATM    7  C12 AOOI    0       1.545  -2.133  -0.358  1.00  0.00           C
HETATM    8  H12 AOOI    0       2.555  -2.439  -0.668  1.00  0.00           H
HETATM    9  H13 AOOI    0       0.875  -2.837  -0.869  1.00  0.00           H
HETATM   10   C7 AOOI    0       1.236  -0.708  -0.875  1.00  0.00           C
HETATM   11  H11 AOOI    0       1.338  -0.813  -1.962  1.00  0.00           H
HETATM   12   C2 AOOI    0      -0.237  -0.338  -0.694  1.00  0.00           C
HETATM   13   C6 AOOI    0      -1.102  -0.728  -1.794  1.00  0.00           C
HETATM   14   O3 AOOI    0      -0.749  -1.333  -2.798  1.00  0.00           O
HETATM   15   C3 AOOI    0      -0.781   0.308   0.394  1.00  0.00           C
HETATM   16   O1 AOOI    0      -0.092   0.661   1.485  1.00  0.00           O
HETATM   17   H1 AOOI    0       0.699   0.070   1.581  1.00  0.00           H
HETATM   18   C4 AOOI    0      -2.188   0.668   0.417  1.00  0.00           C
HETATM   19  C11 AOOI    0      -2.797   1.349   1.488  1.00  0.00           C
HETATM   20   H5 AOOI    0      -2.188   1.643   2.336  1.00  0.00           H
HETATM   21  C10 AOOI    0      -4.156   1.634   1.455  1.00  0.00           C
HETATM   22   H4 AOOI    0      -4.620   2.160   2.284  1.00  0.00           H
HETATM   23   C9 AOOI    0      -4.932   1.233   0.355  1.00  0.00           C
HETATM   24   H3 AOOI    0      -5.997   1.448   0.335  1.00  0.00           H
HETATM   25   C1 AOOI    0      -4.352   0.554  -0.711  1.00  0.00           C
HETATM   26   H2 AOOI    0      -4.932   0.230  -1.569  1.00  0.00           H
HETATM   27   C5 AOOI    0      -2.983   0.280  -0.671  1.00  0.00           C
HETATM   28   O2 AOOI    0      -2.447  -0.399  -1.730  1.00  0.00           O
HETATM   29   C8 AOOI    0       2.237   0.391  -0.488  1.00  0.00           C
HETATM   30  C15 AOOI    0       1.995   1.712  -0.903  1.00  0.00           C
HETATM   31   H6 AOOI    0       1.068   1.948  -1.419  1.00  0.00           H
HETATM   32  C19 AOOI    0       3.447   0.123   0.166  1.00  0.00           C
HETATM   33  H10 AOOI    0       3.681  -0.880   0.505  1.00  0.00           H
HETATM   34  C18 AOOI    0       4.380   1.137   0.405  1.00  0.00           C
HETATM   35   H9 AOOI    0       5.309   0.898   0.917  1.00  0.00           H
HETATM   36  C17 AOOI    0       4.124   2.443  -0.013  1.00  0.00           C
HETATM   37   H8 AOOI    0       4.850   3.232   0.172  1.00  0.00           H
HETATM   38  C16 AOOI    0       2.922   2.726  -0.667  1.00  0.00           C
HETATM   39   H7 AOOI    0       2.706   3.739  -1.000  1.00  0.00           H
CONECT    1    2
CONECT    2    1    3    4    5
CONECT    3    2
CONECT    4    2
CONECT    5    2    6    7
CONECT    6    5
CONECT    7    5    8    9   10
CONECT    8    7
CONECT    9    7
CONECT   10    7   11   12   29
CONECT   11   10
CONECT   12   10   13   15
CONECT   13   12   14   28
CONECT   14   13
CONECT   15   12   16   18
CONECT   16   15   17
CONECT   17   16
CONECT   18   15   19   27
CONECT   19   18   20   21
CONECT   20   19
CONECT   21   19   22   23
CONECT   22   21
CONECT   23   21   24   25
CONECT   24   23
CONECT   25   23   26   27
CONECT   26   25
CONECT   27   18   25   28
CONECT   28   13   27
CONECT   29   10   30   32
CONECT   30   29   31   38
CONECT   31   30
CONECT   32   29   33   34
CONECT   33   32
CONECT   34   32   35   36
CONECT   35   34
CONECT   36   34   37   38
CONECT   37   36
CONECT   38   30   36   39
CONECT   39   38
END''',
    'MZM': '''COMPND    MZM 
AUTHOR    GENERATED BY OPEN BABEL 2.3.90
HETATM    1  N   UNL     1      50.656  41.062  91.081  1.00  0.00           N  
HETATM    2  S   UNL     1      52.283  41.131  90.895  1.00  0.00           S  
HETATM    3  O   UNL     1      52.893  41.831  91.982  1.00  0.00           O  
HETATM    4  O   UNL     1      52.865  39.834  90.926  1.00  0.00           O  
HETATM    5  C   UNL     1      52.669  41.983  89.411  1.00  0.00           C  
HETATM    6  S   UNL     1      53.358  43.339  89.326  1.00  0.00           S  
HETATM    7  C   UNL     1      53.456  43.700  87.843  1.00  0.00           C  
HETATM    8  N   UNL     1      52.884  42.625  87.127  1.00  0.00           N  
HETATM    9  C   UNL     1      52.777  42.542  85.721  1.00  0.00           C  
HETATM   10  N   UNL     1      52.379  41.528  88.122  1.00  0.00           N  
HETATM   11  N   UNL     1      54.022  44.867  87.242  1.00  0.00           N  
HETATM   12  C   UNL     1      54.606  45.978  87.928  1.00  0.00           C  
HETATM   13  O   UNL     1      54.669  46.032  89.135  1.00  0.00           O  
HETATM   14  C   UNL     1      55.159  47.148  87.142  1.00  0.00           C  
HETATM   15  H   UNL     1      50.437  40.576  91.927  1.00  0.00           H  
HETATM   16  H   UNL     1      50.286  41.990  91.124  1.00  0.00           H  
HETATM   17  H   UNL     1      52.289  41.596  85.445  1.00  0.00           H  
HETATM   18  H   UNL     1      52.178  43.386  85.349  1.00  0.00           H  
HETATM   19  H   UNL     1      53.781  42.580  85.274  1.00  0.00           H  
HETATM   20  H   UNL     1      55.553  47.904  87.837  1.00  0.00           H  
HETATM   21  H   UNL     1      55.968  46.798  86.484  1.00  0.00           H  
HETATM   22  H   UNL     1      54.358  47.592  86.533  1.00  0.00           H  
CONECT    1   15   16    2                                            
CONECT    2    4    1    3    5                                       
CONECT    3    2                                                      
CONECT    4    2                                                      
CONECT    5    2    6   10                                            
CONECT    6    5    7                                                 
CONECT    7    6    8   11                                            
CONECT    8    7   10    9                                            
CONECT    9    8   17   18   19                                       
CONECT   10    8    5                                                 
CONECT   11    7   12                                                 
CONECT   12   11   13   14                                            
CONECT   13   12                                                      
CONECT   14   12   20   21   22                                       
CONECT   15    1                                                      
CONECT   16    1                                                      
CONECT   17    9                                                      
CONECT   18    9                                                      
CONECT   19    9                                                      
CONECT   20   14                                                      
CONECT   21   14                                                      
CONECT   22   14                                                      
MASTER        0    0    0    0    0    0    0    0   22    0   22    0
END''',
    'methylazide': '''HEADER    UNCLASSIFIED                            04-Apr-16
TITLE     ALL ATOM STRUCTURE FOR MOLECULE UNK                                   
AUTHOR    GROMOS AUTOMATIC TOPOLOGY BUILDER REVISION 2016-03-31 14:08:37
AUTHOR   2  http://compbio.biosci.uq.edu.au/atb
HETATM    1   N3 _JR3    0      -1.670  -0.314   0.022  1.00  0.00           N
HETATM    2   N2 _JR3    0      -0.592   0.068   0.000  1.00  0.00           N
HETATM    3   N1 _JR3    0       0.512   0.614  -0.024  1.00  0.00           N
HETATM    4   C1 _JR3    0       1.686  -0.290  -0.003  1.00  0.00           C
HETATM    5   H1 _JR3    0       2.566   0.353   0.015  1.00  0.00           H
HETATM    6   H2 _JR3    0       1.683  -0.925   0.889  1.00  0.00           H
HETATM    7   H3 _JR3    0       1.717  -0.919  -0.900  1.00  0.00           H
CONECT    1    2
CONECT    2    1    3
CONECT    3    2    4
CONECT    4    3    5    6    7
CONECT    5    4
CONECT    6    4
CONECT    7    4
END''',
}

OPTIONS = {
    'warfarin': {'total_number_hydrogens': 16, 'net_charge': 0},
}

if __name__ == '__main__':
    for (molecule_name, pdb_str) in PDBS.items():
        molecule = molecule_from_pdb_str(pdb_str, name=molecule_name)
        if molecule.name in {'warfarin'}:
            print(molecule.get_all_tautomers(**OPTIONS[molecule_name] if molecule_name in OPTIONS else {}))
        else:
            print(molecule.assign_bond_orders_and_charges_with_ILP(enforce_octet_rule=True))
        print(molecule.write_graph(molecule_name, output_size=(int(2100 / 1.5), int(2970 / 1.5))))
        if molecule_name == 'warfarin':
            print(molecule)
        print()

    molecule = Molecule([Atom(index=1, element='C', valence=3, capped=True, coordinates=None), Atom(index=2, element='C', valence=3, capped=True, coordinates=None)], [(1,2)])
    print(molecule.get_all_tautomers())
    print(molecule.write_graph('ethene', output_size=(200, 200)))

    molecule = Molecule([Atom(index=1, element='H', valence=1, capped=True, coordinates=None), Atom(index=2, element='O', valence=1, capped=True, coordinates=None)], [(1, 2)], netcharge=0, name='hydroxyl_radical')
    print(molecule.assign_bond_orders_and_charges_with_ILP())
    print(molecule.write_graph('', output_size=(200, 200)))
