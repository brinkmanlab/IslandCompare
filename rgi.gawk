BEGIN { FS=OFS="\t"; print "##gff-version 3" }
NR == 1 { split( $0, tags ); next }
{ print $2,"RGI-CARD","gene",$3,$4,$8,$5,".","Name=AMR Gene;Alias=ARO:"$11";"tags[6]"="$6";"tags[7]"="$7";"tags[9]"="$9";"tags[10]"="$10";"tags[12]"="$12";"tags[13]"="$13";"tags[14]"="$14";"tags[16]"="$16";"tags[17]"="$17";"tags[15]"="$15";" }
