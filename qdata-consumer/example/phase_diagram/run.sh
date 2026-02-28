../../build/sample phasediagram.ini
../../build/learn phasediagram.clone.h5 --noInfT
../../build/segregate-phases phasediagram.out.h5 --weight=lorentzian
python3 plot_phasediagram.py
