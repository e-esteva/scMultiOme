getCellDex=function(obj,species,assay='RNA'){
  
  trained.sce=SingleCellExperiment(list(counts=obj@assays[[assay]]@data))
  
  if(species == 'mouse'){
    ref <- ImmGenData()
    
    
    broad_anno <- SingleR(test = trained.sce, ref = ref, assay.type.test=1,
                          labels = ref$label.main)
    write.csv(table(broad_anno$labels),'broad-annotations-summary.csv')
    
    obj$singleR_ImmgenLabels_broad=broad_anno$labels
    Idents(obj)=obj$singleR_ImmgenLabels_broad
    DefaultAssay(obj)=assay
    p=DimPlot(obj,label = T)
    ggsave('ImmGenAnnotations-broad.pdf',p)
    DefaultAssay(obj)='RNA'
    obj = obj %>% NormalizeData(obj)
    broad_anno.markers=FindAllMarkers(obj,only.pos=T)
    write.csv(broad_anno.markers,'broad-anno-markers.csv')
    top10_markers=broad_anno.markers %>% data.frame() %>% group_by(cluster) %>% top_n(10,avg_log2FC)
   
    DefaultAssay(obj)=assay

    if(dim( obj@assays[[assay]]@scale.data )[1] ==0){

        obj=ScaleData(obj)
    }    
    
    h=DoHeatmap(obj,features=top10_markers$gene)
    ggsave('broad-annotations-top10-markers-Heatmap.pdf',h)
    
    
    fine_anno <- SingleR(test = trained.sce, ref = ref, assay.type.test=1,
                         labels = ref$label.fine)
    write.csv(table(fine_anno$labels),'fine-annotations-summary.csv')
    
    obj$singleR_ImmgenLabels_fine=fine_anno$labels

    saveRDS(obj,'broad-fine-annos-obj.rds')

    Idents(obj)=obj$singleR_ImmgenLabels_fine
    DefaultAssay(obj)=assay
    p=DimPlot(obj,label = T)
    ggsave('ImmGenAnnotations-fine.pdf',p,height=20,width=20)
    
    DefaultAssay(obj)='RNA'
    fine_anno.markers=FindAllMarkers(obj,only.pos=T)
    write.csv(fine_anno.markers,'fine-anno-markers.csv')
    DefaultAssay(obj)=assay
    
    top10_markers=fine_anno.markers %>% data.frame() %>% group_by(cluster)	%>% top_n(10,avg_log2FC)
    ranges_	= seq(1,length(unique(top10_markers$cluster)),5)
    if(length(ranges_) > 1){
      if(max(ranges_) < length(unique(top10_markers$cluster))){
        ranges_[length(ranges_)]=length(unique(top10_markers$cluster))
      }

      ranges_[1]=0

      for(i in seq(2,length(ranges_))){

	cluster.tmp=unique(top10_markers$cluster)[seq(ranges_[i-1]+1,ranges_[i])]

        tmp=subset(obj,id=cluster.tmp)

	features.tmp=top10_markers$gene[top10_markers$cluster%in%cluster.tmp]
	
	if(length(features.tmp) > 15){
                h=DoHeatmap(tmp,features=features.tmp)+theme(text=element_text(size=4))
          }else{
                h=DoHeatmap(tmp,features=features.tmp)
          }

        ggsave(glue('fine-annotations-top10-markers-Heatmap-{i-1}.pdf'),h)
      }
    }else{
      h=DoHeatmap(obj,features=top10_markers$gene)
      ggsave('fine-annotations-top10-markers-Heatmap.pdf',h)
      
    }
  }else{
    
    ref1 <- HumanPrimaryCellAtlasData()
    ref2 <-  DatabaseImmuneCellExpressionData()
    
    refs=list(ref1,ref2)
    ref_names=c('HumanPrimaryCellAtlas','DatabaseImmuneCellExpression')
    
    for(r in seq(length(refs))){
      
      ref = refs[[r]]
      ref.dir.tmp=ref_names[r]
      dir.create(ref.dir.tmp)

      broad_anno <- SingleR(test = trained.sce, ref = ref, assay.type.test=1,
                            labels = ref$label.main)
      
      write.csv(table(broad_anno$labels),glue('{ref.dir.tmp}/broad-annotations-summary.csv'))
      
      obj[[glue('{ref.dir.tmp}_broad')]]=broad_anno$labels
      Idents(obj)=obj[[glue('{ref.dir.tmp}_broad')]]
      DefaultAssay(obj)='RNA'
      p=DimPlot(obj,label = T)
      ggsave(glue('{ref.dir.tmp}/annotations-broad.pdf'),p)
      broad_anno.markers=FindAllMarkers(obj,only.pos=T)
      
      write.csv(broad_anno.markers,glue('{ref.dir.tmp}/broad-anno-markers.csv'))
      
      top10_markers=broad_anno.markers %>% data.frame() %>% group_by(cluster) %>% top_n(10,avg_log2FC)
      obj=ScaleData(obj)
      
      h=DoHeatmap(obj,features=top10_markers$gene)
      ggsave(glue('{ref.dir.tmp}/broad-annotations-top10-markers-Heatmap.pdf'),h)
      
      
      fine_anno <- SingleR(test = trained.sce, ref = ref, assay.type.test=1,
                           labels = ref$label.fine)
      write.csv(table(fine_anno$labels),glue('{ref.dir.tmp}/fine-annotations-summary.csv'))
      
      obj[[glue('{ref.dir.tmp}_fine')]]=fine_anno$labels
      saveRDS(obj,'broad-fine-annos-obj.rds')
      Idents(obj)=obj[[glue('{ref.dir.tmp}_fine')]]
      DefaultAssay(obj)='RNA'
      p=DimPlot(obj,label = T)
      ggsave(glue('{ref.dir.tmp}/annotations-fine.pdf'),p,width=20,height=20)
      
      fine_anno.markers=FindAllMarkers(obj,only.pos=T)
      write.csv(fine_anno.markers,glue('{ref.dir.tmp}/fine-anno-markers.csv'))
      
      top10_markers=fine_anno.markers %>% data.frame() %>% group_by(cluster)	%>% top_n(10,avg_log2FC)
      ranges_	= seq(1,length(unique(top10_markers$cluster)),5)
      if(length(ranges_) > 1){
        if(max(ranges_) < length(unique(top10_markers$cluster))){
          ranges_[length(ranges_)]=length(unique(top10_markers$cluster))
        }
        ranges_[1]=0
        for(i in seq(2,length(ranges_))){
	
	  cluster.tmp=unique(top10_markers$cluster)[seq(ranges_[i-1]+1,ranges_[i])]

          tmp=subset(obj,id=cluster.tmp)
	 
          features.tmp=top10_markers$gene[top10_markers$cluster%in%cluster.tmp]

	  if(length(features.tmp) > 15){
		h=DoHeatmap(tmp,features=features.tmp)+theme(text=element_text(size=4))
	  }else{
		h=DoHeatmap(tmp,features=features.tmp)
	  }


          ggsave(glue('{ref.dir.tmp}/fine-annotations-top10-markers-Heatmap-{i-1}-.pdf'),h)
        }
      }else{
        h=DoHeatmap(obj,features=top10_markers$gene)
        ggsave('{ref.dir.tmp}/fine-annotations-top10-markers-Heatmap.pdf',h)
        
      }
      
    }
    
  }
  
  
}
