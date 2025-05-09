##rbhXpress.sh -a $work_dir/db/$species[$i].prot.fa -b $work_dir/db/$species[$j].prot.fa -o $species[$i].$species[$j].tab -t 8

library(macrosyntR)
library(ggplot2)
species<-c("Alref","Bb","Bf","Bj","Bl","Ea","Gg","Hs","Lo","Lv","Mm","Oa","Pf","Py","Xt")

for(i in 1:14){
	a<-i+1
	for(j in a:15){  
		print(paste("processing",species[i],species[j],sep=" "))
	  file1=paste("F:/chrom/",species[i],".",species[j],".tab",sep="")
	  file2=paste("F:/chrom/",species[i],".bed",sep="")
	  file3=paste("F:/chrom/",species[j],".bed",sep="")
	  my_orthologs <- load_orthologs(orthologs_table = file1,
                                   sp1_bed = file2,
                                   sp2_bed = file3)
                                   
    p2 <- plot_oxford_grid(my_orthologs,
                           sp1_label = species[i],
                           sp2_label = species[j],
                           reorder = TRUE,
                           color_by = "clust")
                           
    my_orthologs_reordered <- reorder_macrosynteny(my_orthologs)
    my_macrosynteny <- compute_macrosynteny(my_orthologs_reordered)
    p3 <- plot_macrosynteny(my_macrosynteny)                      

    output1=paste("F:/dotplot/",species[i],".",species[j],".dotplot.pdf",sep="")
    output2=paste("F:/dotplot/",species[i],".",species[j],".test.pdf",sep="")
    output3=paste("F:/dotplot/",species[i],".",species[j],".test.tbl",sep="")
    
    ggsave(p2,file=output1)
    ggsave(p3,file=output2)
    
    write.table(my_macrosynteny,output3)
	}
}
