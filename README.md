# Comandos para extrair resultados de uma simulação utilizando o programa GROMACS

## Após termino da dinâmica:

Texto adaptado do [Manual do GROMACS](https://www.gromacs.org/Documentation_of_outdated_versions/Terminology/Periodic_Boundary_Conditions)

Condições de contorno periódicas (do inglês, *Periodic Boundary Condition*, ou PBC) são utilizadas em simulações de dinâmica molecular para evitar problemas com efeitos de contorno causados por tamanho finito e tornar o sistema mais infinito, às custas de possíveis efeitos de periodicidade.
Iniciantes visualizando uma trajetória às vezes pensam que estão observando um problema quando:
1. A(s) molécula(s) não ficam no centro da caixa;
2. Parece que partes da molécula(s) se difundem fora da caixa;
3. Buracos são criados;
4. Moléculas quebradas aparecem;
5. Sua caixa de simulação era um dodecaedro rômbico ou octaedro cúbico, mas parece um cubo inclinado após a simulação;
6. Ligações malucas em toda a célula de simulação aparecem.

Este não é um problema ou erro que está ocorrendo, é o que você deve esperar.
A existência de PBC significa que qualquer átomo que sai de uma caixa de simulação, digamos, pela face direita, entra na caixa de simulação pela face esquerda. No exemplo de uma proteína grande, se você olhar para a face da caixa de simulação que é oposta àquela de onde a proteína está saindo, então um buraco no solvente será visível. A razão pela qual as moléculas se movem de onde estavam inicialmente localizadas dentro da caixa é (para a grande maioria das simulações) que elas podem se difundir livremente. E é o que eles fazem, eles não são mantidos em um local mágico da caixa.

Um exemplo animado está a baixo. Perceba que a caixa do meio pode ser considerada a principal, enquanto "cópias" desta caixa são colocadas em volta.

<p align="center">
  <img src="https://i.makeagif.com/media/11-11-2017/f9YXHQ.gif" alt="animated" />
</p>

Para retirar as condições de contorno periódicas, utilize o comando:

```
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -ur compact 
```
* Selecionar System (0)

OBS.: Caso esteja utilizando mais de uma estrutura no seu sistema, considere utilizar `nojump` como argumento na **flag** `-pbc`


## Root-Mean-Square Deviation (RMSD)

Texto retirado do [Portal BIOINFO](https://bioinfo.com.br/dinamica-molecular-como-mostrar-um-filme-completo-em-uma-folha-de-papel/):

>Ao avaliar um sistema que está em movimento, você precisa de um gráfico que represente sua mobilidade. Uma forma de fazer isso é utilizar o desvio quadrático-médio das distâncias dos átomos (ou do inglês *Root-Mean-Aquare Deviation* ou somente RMSD). Nesse tipo de gráfico é feita uma comparação frame a frame da variação das distâncias. Quando o gráfico alcança o platô, ou seja, não tiver mais tantas variações, pode-se dizer que o sistema entrou em equilíbrio, indicando que a proteína, por exemplo, não apresenta mais tantas modificações estruturais.

```
gmx rms -s md.tpr -f md_noPBC.xtc -o rmsd.xvg -tu ns
```
* Selecionar Backbone (4) e Backbone (4) ou um index desejado.

## Raio de Giro

Texto retirado do [Portal BIOINFO](https://bioinfo.com.br/dinamica-molecular-como-mostrar-um-filme-completo-em-uma-folha-de-papel/):

>Bastante usado para estudo de enovelamento de peptídeos e proteínas. O raio de giro está relacionado ao deslocamento do centro de massa da proteína em relação a um eixo. Simplificando, quanto mais volumosa uma proteína, maior o seu raio de giro. Nesse caso, podemos dizer que quando uma proteína está desestruturando, o seu raio de gira aumenta.

```
gmx gyrate -s md.tpr -f md_noPBC.xtc -o gyrate.xvg
```

* Selecionar Protein (1) ou um index desejado.

## Root-Mean-Square Fluctuation (RMSF)

A flutuação quadrática média (*do ingês, Root-Mean-Square Fluctuation*, ou RMSF) de uma estrutura é a média de tempo do RMSD. O RMSD quantifica o quanto uma estrutura diverge de uma referência ao longo do tempo, o RMSF revela quais áreas do sistema são mais móveis.

```
gmx rmsf -s md.tpr -f md_noPBC.xtc -o rmsf.xvg -res (Backbone)
```
* Selecionar Backbone (4) ou um index desejado.

## Para recomeçar a dinâmica

```
gmx convert-tpr -s md.tpr -extend 50000 -o md_150.tpr 
```
* Estender o tempo de simulação por 50ns. Perceba que o tempo está em ps (50000)

GROMACS utiliza o arquivo `md.cpt` para conter todas as informações necessárias para continuar a simulação:

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

Caso tenha optado por separar os arquivos da dinâmica, você pode junta-los utilizando:

```
gmx trjcat -f md_1.xtc md_2.xtc -o trajectory.xtc
```
```
gmx trjcat -f md_1.trr md_2.trr -o trajectory.trr
```

## Utilizar Algoritmos de clusterização de estruturas:

O GROMACS possui uma ferramenta para poder agrupar estruturas utilizando diferentes métodos. As informações abaixo foram retiradas do [Manual](https://manual.gromacs.org/documentation/2021/onlinehelp/gmx-cluster.html)

* `linkage`: adiciona uma estrutura a um cluster quando sua distância a qualquer elemento do cluster for menor que um valor de corte.

* `jarvis-patrick`: adiciona uma estrutura a um cluster quando esta estrutura é vizinha a outra estrutura no cluster tiverem pelo menos P vizinhos em comum. Os vizinhos de uma estrutura são as M estruturas mais próximas ou todas as estruturas dentro de um valor de corte.

* `monte-carlo`: o algoritmo reordena a matriz de RMSD de forma que a ordem dos quadros use os menores incrementos possíveis. Com isso, é possível fazer uma animação suave indo de uma estrutura para outra com o maior RMSD possível (por exemplo) entre elas, porém as etapas intermediárias devem ser as menores possíveis. As aplicações podem ser para visualizar um potencial de conjunto de força média de simulações ou uma simulação de tração. Obviamente, o usuário deve preparar bem a trajetória (por exemplo, não sobrepondo quadros). O resultado final pode ser inspecionado visualmente olhando para o arquivo `.xpm` da matriz, que deve variar suavemente de baixo para cima.

* `diagonalization`: utiliza a diagonal da matriz de RMSD.

* `gromos`: usa um algoritmo descrito em Daura et al. (Angew. Chem. Int. Ed. 1999, 38, pp 236-240). Conta o número de vizinhos utilizando um valor de corte, pega a estrutura com maior número de vizinhos com todos os seus vizinhos como um cluster e a retira do `pool` de clusters. O algoritmo repete esta etapa para todas as estruturas restantes.

Arquivos necessários: `md_150.tpr`, `md.tpr`, `md.trr`, `md_noPBC.xtc`, `topol.top`

Para fazer a clusterização:

```
gmx cluster -f md.trr -s md.tpr -cutoff 0.1 -dist -sz -ntr -clid -wcl 0 -method linkage
```
OBS.: Selecione Backbone (4) e Protein (1) para utilizar os átomos principais como entrada para os cálculos e o programa retornar arquivos `.pdb` com todos os átomos.

## Para criar index customizáveis

Adaptado do [Manual](https://manual.gromacs.org/archive/5.0.2/programs/gmx-make_ndx.html):

>Os grupos de índice são necessários para quase todos os programas GROMACS. Todos esses programas podem gerar grupos de índice padrão. Você SÓ precisa usar `gmx make_ndx` quando precisar de grupos de índice ESPECIAIS. Há um grupo de índice padrão para todo o sistema, 9 grupos de índice padrão para proteínas e um grupo de índice padrão é gerado para cada outro nome de resíduo. Quando nenhum arquivo de índice é fornecido, o `gmx make_ndx` irá gerar os grupos padrão. Com o editor de índice, você pode selecionar nomes e números de átomos, resíduos e cadeias. Quando um arquivo de entrada de execução é fornecido, você também pode selecionar o tipo de átomo. Você pode usar NOT, AND e OR, você pode dividir grupos em cadeias, resíduos ou átomos. Você pode excluir e renomear grupos.

```
gmx make_ndx -f md.tpr -o n.ndx
```

OBS.: Você pode combinar opções, como por exemplo, selecionar somente os átomos do *backbone* de resíduos específicos, digitando `r 59-171 & 4`.

## Criar PDB com valores de B-Factor:

É possível renderizar cores em uma proteína a partir dos valores de RMSF, o que permite identificação visual das regiões de alta flutuação. Isso é comumente feito configurando os valores do fator de temperatura (também conhecido como b-factor), o que pode ser visualizado utilizando programas como VMD ou o UCSF Chimera.

```
gmx rmsf -s md.tpr -f md_noPBC.xtc -oq rmsf.pdb -res
```

* Selecionar Backbone (4), caso deseje somente os átomos principais ou Protein (1) para poder ter todos os átomos e estruturas secundárias.

## Para medir a distancia entre 2 grupos:

```
gmx pairdist -f md_noPBC.xtc -s md.tpr -o paidist.xvg -tu ns
```

## DSSP:

Estre programa lê um arquivo de trajetória e calcula a estrutura secundária para cada período de chamada do programa `DSSP`. Se você não tiver o programa `DSSP`, acesse <https://swift.cmbi.umcn.nl/gv/dssp/> ou utilize o disponível neste repositório.

Primeiramente é necessário colocar o arquivo `dssp` dentro de **/usr/local/bin** (ou em alguma pasta no seu $PATH) e em seguida transformar em executável: `chmod +x /usr/local/bin/dssp`

```
gmx do_dssp -f trajetória -s *.tpr -o dssp.xpm -tu ns -sc scount.xvg -sss HE
```

```
gmx xpm2ps -f dssp.xpm -o dssp.eps -di in.m2p 
```

* Para alterar os parametros do gráfico, editar os arquivos `.xpm` e `.m2p`.

## Solvent-Accessible Surface Area (SASA):

O `gmx sasa` calcula as áreas acessível pelo solvente utilizando o algoritmo de Eisenhaber F, Lijnzaad P, Argos P, Sander C, & Scharf M (1995) J. Comput. Chem. 16, 273-284

```
gmx sasa -f md_noPBC.xtc -s md.tpr -n n.ndx -o sasa.xvg -or sasa_perresidue.xvg -tu ns -pbc yes
```

OBS.: Lembrar de criar um arquivo index caso seja preciso selecionar resíduos específicos

## HBonds

`gmx hbond` calcula e analisa ligações de hidrogênio. As ligações de hidrogênio são determinadas com base nos pontos de corte para o ângulo `Hidrogênio - Doador - Receptor`  e a distância `Doador - Receptor` (ou Hidrogênio - Receptor usando `-noda`). Os grupos `OH` e `NH` são considerados **doadores**, `O` é sempre um **receptor**, N é um **receptor** por padrão.

```
gmx hbond -f md_noPBC.xtc -s md.tpr -tu ns
```

## Editar aquivos .xvg

Esta série de comandos ira converter o arquivo `.xvg` gerado pelo GROMACS em arquivos `.txt` que pode ser utilizado para fazer graficos em R ou Python.

* Raio de Giro:

```
cat gyrate.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > gyrate.txt && sed -i 's/@TYPE xy/time rg rx ry rz/' gyrate.txt
```

* SASA:

```
cat sasa.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > sasa.txt && sed -i 's/@TYPE xy/time value/' sasa.txt
```


* SASA por resíduo:

```
cat sasa_perresidue.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > sasa_perresidue.txt && sed -i 's/@TYPE xy/time avg sd/' sasa_perresidue.txt
```

* RMSD:

```
cat rmsd.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > rmsd.txt && sed -i 's/@TYPE xy/time value/' rmsd.txt
``` 

* RMSF:

``` 
cat rmsf.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > rmsf.txt && sed -i 's/@TYPE xy/time value/' rmsf.txt
```

* Pairdist:

```
cat pairdist.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > pairdist.txt && sed -i 's/@TYPE xy/time value/' pairdist.txt 
```

* HBonds:

```
cat hbnum.xvg | awk '{if ($1 != "@" && $1 != "#") print}' > hbnum.txt && sed -i 's/@TYPE xy/time total belowcutoff/' hbnum.txt 
```
