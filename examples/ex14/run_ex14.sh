#!/bin/bash

EXEMDFF=mdff.x

mdff=true

if $mdff; then
	echo "running mdff test (ex14)"
	rm -rf mdff 
	mkdir -p mdff
	cd mdff
		cp ../config/control_* .
		echo "NVE velocities rescaled lennard-jones"
		cp ../config/POSFF.Ar POSFF 
		$EXEMDFF control_md_full_run_lj.F > stdout_md_full_run_lj
		echo "full run finish"
		mv OSZIFF OSZIFF.full_lj
		$EXEMDFF control_md_splited_run1_lj.F > stdout_md_splited_run1_lj
		echo "split 1 finish"
		mv OSZIFF OSZIFF.run1_lj
		$EXEMDFF control_md_splited_run2_lj.F > stdout_md_splited_run2_lj
		echo "split 2 finish"
		mv OSZIFF OSZIFF.run2_lj
		echo "NVE velocities rescaled PIM"
		cp ../config/POSFF.SiO2 POSFF
		$EXEMDFF control_md_full_run_pim.F > stdout_md_full_run_pim
		echo "full run finish"
		mv OSZIFF OSZIFF.full_pim
		$EXEMDFF control_md_splited_run1_pim.F > stdout_md_splited_run1_pim
		echo "split 1 finish"
		mv OSZIFF OSZIFF.run1_pim
		$EXEMDFF control_md_splited_run2_pim.F > stdout_md_splited_run2_pim
		echo "split 2 finish"
		mv OSZIFF OSZIFF.run2_pim
		for c in full_lj run1_lj run2_lj full_pim run1_pim run2_pim
		do
			poszi -n -i OSZIFF.$c
			mv etot_l etot_l.$c
			mv temp_l temp_l.$c
		done
	cd ..
fi

cat > plot.lj << eof
#!/usr/bin/gnuplot -persist
reset
set xlabel "# step"
set ylabel "E_pot (eV)"
p 'mdff/etot_l.full_lj' u 3:15 w l  lc 3 title "FULL run",'mdff/etot_l.run1_lj' u 3:15  w p lc 1 title "RUN 1",'mdff/etot_l.run2_lj' u 3:15 w p lc 2  title "RUN 2" 
eof
chmod u+x plot.lj
./plot.lj

cat > plot.pim << eof
#!/usr/bin/gnuplot -persist
reset
set xlabel "# step"
set ylabel "E_pot (eV)"
p 'mdff/etot_l.full_pim' u 3:15 w l  lc 3 title "FULL run",'mdff/etot_l.run1_pim' u 3:15  w p lc 1 title "RUN 1",'mdff/etot_l.run2_pim' u 3:15 w p lc 2  title "RUN 2" 
eof
chmod u+x plot.pim
./plot.pim

