#/bin/bash -i
#
car_list=" Pt_surf_100_clean.car Pt_surf_100_H2_only.car Pt_surf_100_hollow.car Pt_surf_100_top.car Pt_surf_100_top_opt.car Pt_surf_111_clean.car Pt_surf_111_H2_only.car Pt_surf_111_hollow.car Pt_surf_111_top.car "

for car_file in $car_list; do

   if [ -e "$car_file" ]; then

      build_poscar="0"
      stem=`echo $car_file | cut -f1 -d.`

      if [ -e "$stem" ]; then
         echo Directory for $stem already exists...
         if [ -e "$stem"/POSCAR_latest ]; then
           echo There is already a POSCAR_latest in this directory so doing nothing..
         else
           echo Just making the POSCAR_latest file..
           build_poscar="1"
         fi
      else
         mkdir $stem
         build_poscar="1"
      fi
#
# Build the POSCAR file, if required
#
      echo build_poscar = $build_poscar
      if [ $build_poscar = "1" ]; then
         echo Building POSCAR_$stem
         echo master $car_file > inter.inp
         echo zsort >> inter.inp
         echo need poscar POSCAR_$stem >> inter.inp
         inter_vasp inter.inp > "$stem"_inter.out
         cp POSCAR_$stem $stem/POSCAR_latest
      else
         echo NOT building POSCAR_$stem
      fi
#
      cat hawk_opt_blank.sh | sed s/BLANK/$stem/ > $stem.sh

    else
      echo $car_file does not exist in this directory, check list in script
  fi

done

