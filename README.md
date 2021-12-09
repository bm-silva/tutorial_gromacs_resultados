# Comandos para extrair resultados de uma simulação utilizando o programa GROMACS

## Após termino da dinâmica:

```
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -ur compact 
```
* Selecionar System (0)

OBS.: Caso esteja utilizanod mais de uma estrutura no seu sistema, considere utilizar `nojump` como argumento na **flag** `-pbc`


## Root-Mean-Square Deviation (RMSD)

```
gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns
```
* Selecionar Backbone (4) e Backbone (4) ou um index desejado.

## Raio de Giro

```
gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg
```

* Selecionar Protein (1) ou um index desejado.

## Root-Mean-Square Fluctuation (RMSF)

```
gmx rmsf -s md.tpr -f md_noPBC.xtc -o rmsf.xvg -res (Backbone)
```
* Selecionar Backbone (4) ou um index desejado.

## Para recomeçar a dinâmica


```
gmx convert-tpr -s md.tpr -extend 50000 -o md_150.tpr 
```
* Estender o tempo de simulação por 50ns. Perceba que o tempo está em ps (50000)

1. Caso a sua dinâmica tenha parado antes do tempo, você pode reinicia-la digitando:
```
gmx mdrun -s md.tpr -cpi md.cpt -v -deffnm md 
```
2. Para recomeçar a dinâmica e salvar todos os novos arquivos em arquivos separados:
```
gmx mdrun -s md_150.tpr -cpi md.cpt -noappend -v -deffnm md_150 
```
3.  Para recomeçar a dinâmica e salvar todos os novos arquivos dentro dos arquivos ja existentes:
```
gmx mdrun -s md_150.tpr -cpi md.cpt -v -append -deffnm md 
```

## Para juntar arquivos de trajetoria:

```
gmx trjcat -f md_1.xtc md_2.xtc -o trajectory.xtc
```
```
gmx trjcat -f md_1.trr md_2.trr -o trajectory.trr
```

## Utilizar Algoritmos de clusterização de estruturas:

Arquivos necessários: `md_150.tpr`, `md.tpr`, `md.trr`, `md_noPBC.xtc`, `topol.top`

* Para fazer a clusterização:

```
gmx cluster -f md.trr -s md.tpr -cutoff 0.1 -dist -sz -ntr -clid -wcl 0 -method linkage (usa backbone e protein)
```
OBS.: Selecione Backbone (4) e Protein (1) para utilizar os átomos principais como entrada para os cálculos e o programa retornar arquivos `.pdb` com todos os átomos.

## Para criar index customizáveis

```
gmx make_ndx -f md.tpr -o n.ndx
```

OBS.:Você pode combinar opções, como por exemplo, selecionar somente os átomos do **backbone** de resíduos específicos, digitando `r 59-171 & 4`.

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
