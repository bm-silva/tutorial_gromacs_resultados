# Comandos para extrair resultados de uma simulação utilizando o programa GROMACS

Após termino da dinâmica:

	gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -ur compact (selecionar System)

	-pbc nojump para mais de uma chain

	gmx trjconv -f md_noPBC.xtc -s md.tpr -o md_noPBCfit.xtc -fit rot+trans (selecionar 1 e 1:usar para fazer análises envolvendo somente a proteína)


RMSD: gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns (Backbone e Backbone)

Gyrate: gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg (Protein)

RMSF: gmx rmsf -s md.tpr -f md_noPBC.xtc -o rmsf.xvg -res (Backbone)

Para recomeçar a dinâmica

	gmx convert-tpr -s md.tpr -extend 50000 -o md_150.tpr (extende por 50ns, tempo em ps)

	gmx mdrun -v -deffnm md_150 -cpi md.cpt -noappend (para fazer em arquivos separados)

	gmx mdrun -s md_150.tpr -cpi md.cpt -v -append -deffnm md (para fazer tudo junto)

Para juntar arquivos de trajetoria:

	gmx trjcat -f md_1.xtc md_2.xtc -o trajectory.xtc

	gmx trjcat -f md_1.trr md_2.trr -o trajectory.trr

Arquivos para cluster:

	md_150.tpr, md.tpr, trajectory.trr, trajectory.xtc, topol.top

Para fazer cluster:

	gmx cluster -f md_150.part0002.trr -s md_150.tpr -b 100000 -e 140000 -cutoff 0.1 -dist -sz -ntr -clid -wcl 0 -method linkage (usa backbone e protein)

Para criar o index e selecionar o rmsf de aminoácidos específicos

	gmx make_ndx -f md_150.tpr -o n.ndx (r resInicial-resFinal)
	gmx make_ndx -f md.tpr -o n.dx (r 59-171 & 4) para selecionar somente backbone

Criar PDB com B-Score usando o index:

	gmx rmsf -s md.tpr -f md_noPBC.xtc -oq rmsf.pdb -res -n n.ndx (seleciona o index escolhido)

Para medir a distancia entre 2 residuos, ou 2 index:

	gmx pairdist -f trajectory_500_noPBC.xtc -s md_500.tpr -n n2.ndx -o teste_rmsd_sidechain.xvg -tu ns

Para utilizar o DSSP:
	-Primeiramente colocar o arquivo dentro de /usr/local/bin
	-Transformar em executável
	
	gmx do_dssp -f trajetória -s *.tpr -o dssp.xpm -tu ns -sc scount.xvg -sss HE
	gmx xpm2ps -f dssp.xpm -o dssp.eps -di in.m2p (Para alterar os parametros do grafico editar os arquivos .xpm e .m2p)

Para realizar análise de SASA:
	
	gmx sasa -f md_noPBC.xtc -s md.tpr -n n.ndx -o sasa.xvg -or sasa_perresidue.xvg -tu ns -pbc yes

	-Lembrar de criar um arquivo index caso seja preciso selecionar resíduos específicos
