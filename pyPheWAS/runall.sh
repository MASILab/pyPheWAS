#/bin/bash
path="/Users/Nabar/Documents/MASI/Spring_2016_MASI/Saha_Cutting_Data/Cutting_20160318/dirs"
groupext="_group.csv"
icdext="_icd9.csv"
featureheader='feature_matrix_'
regheader='regressions_'
plot='_plot.png'
plotwith='_plotimb.png'
for D in $path/*/; do
	if [ -d ${D} ]; then
		# path to use for pyPhewas is ${D}
		base=$(basename $D)
		groupfile=$base$groupext
		icdfile=$base$icdext
		echo ${D}$groupfile
		head -1 ${D}/$groupfile
		fm=$featureheader$groupfile
		sf=$regheader$groupfile
		p=$base$plot
		pw=$base$plotwith
		./pyPhewasLookup --path=${D} --phenotype=$icdfile --group=$groupfile --reg_type='log'
		./pyPhewasModel --path=${D} --feature_matrix=$fm --covariates='genotype' --reg_type='log' --group=$groupfile
		./pyPhewasPlot --path=${D} --statfile=$sf --imbalance=False --thresh_type='bon' --outfile=$p
		./pyPhewasPlot --path=${D} --statfile=$sf --imbalance=True --thresh_type='bon' --outfile=$pw
	fi
done

