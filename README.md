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

OBS.: Você pode combinar opções, como por exemplo, selecionar somente os átomos do *backbone* de resíduos específicos, digitando `r 59-171 & 4`.

## Criar PDB com valores de B-Score:

```
gmx rmsf -s md.tpr -f md_noPBC.xtc -oq rmsf.pdb -res
```

## Para medir a distancia entre 2 grupos:

```
gmx pairdist -f md_noPBC.xtc -s md.tpr -o paidist.xvg -tu ns
```

## DSSP:

Primeiramente é necessário colocar o arquivo `dssp` dentro de **/usr/local/bin** e em seguida transformar em executável: `chmod +x /usr/local/bin/dssp`

```
gmx do_dssp -f trajetória -s *.tpr -o dssp.xpm -tu ns -sc scount.xvg -sss HE
```

```
gmx xpm2ps -f dssp.xpm -o dssp.eps -di in.m2p 
```

* Para alterar os parametros do gráfico, editar os arquivos `.xpm` e `.m2p`.

## Solvent-Accessible Surface Area (SASA):
	
```
gmx sasa -f md_noPBC.xtc -s md.tpr -n n.ndx -o sasa.xvg -or sasa_perresidue.xvg -tu ns -pbc yes
```

OBS.: Lembrar de criar um arquivo index caso seja preciso selecionar resíduos específicos

## Editar aquivos .xvg

Esta série de comandos ira converter o arquivo `.xvg` gerado pelo GROMACS em arquivos `.txt` que pode ser utilizado para fazer graficos em R ou Python.

* Raio de Giro:

```
cat gyrate.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > gyrate.txt && sed -i 's/@TYPE xy/time rg x y z/' gyrate.txt
```

* SASA:

```
cat sasa.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > sasa.txt && sed -i 's/@TYPE xy/x y/' sasa.txt
```


* SASA por resíduo:

```
cat sasa_perresidue.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > sasa_perresidue.txt && sed -i 's/@TYPE xy/x y z/' sasa_perresidue.txt
```

* RMSD:

```
cat rmsd.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > rmsd.txt && sed -i 's/@TYPE xy/x y/' rmsd.txt
``` 

* RMSF:

``` 
cat rmsf.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > rmsf.txt && sed -i 's/@TYPE xy/x y/' rmsf.txt
```

* Pairdist:

```
cat pairdist.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > pairdist.txt && sed -i 's/@TYPE xy/x y/' pairdist.txt 
```

* HBonds:

```
cat hbnum.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > hbnum.txt && sed -i 's/@TYPE xy/x y z/' hbnum.txt 
```
