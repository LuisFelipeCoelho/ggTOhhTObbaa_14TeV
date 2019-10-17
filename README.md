# ggTOhhTObbaa_14TeV

Correr no ficheiro **/ExpertMode/Build**
```
source setup.sh
make
make clean
./MadAnalysis5job ../Input/in.txt
```

Para correr é preciso definir um arquivo de texto com os "paths" dos eventos, 
neste caso o arquivo está em **/ExpertMode/Input/in.txt**

e os eventos estão em:

**/Events_caixa/run_01_decayed_1/unweighted_events.lhe.gz**

**/Events_triangulo/run_01_decayed_1/unweighted_events.lhe.gz**

**/Events_caixa_e_triangulo/run_02_decayed_1/unweighted_events.lhe.gz**
