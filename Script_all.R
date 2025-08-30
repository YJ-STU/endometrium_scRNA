#Metacell analysis
  library(foreach)
  library(pheatmap)
  library(metacell)
  library(tgconfig)
  library(Matrix)
  library(tgstat)
  setwd("/dir/to/mtx/files/")
  set_param("mc_plot_device",'pdf', "metacell")
  set_param("scm_spike_regexp","^ERCC-","metacell")
  set_param("scm_mc_mark_k_per_clust",100,"metacell") #default: 5
  set_param("scm_mc_mark_min_gene_cov",0.3,"metacell") # default: 0.25
  set_param("scm_mc_mark_min_gene_fold",2,"metacell") # default: 1.5
  set_param("mc_cores",10,"metacell")
  
  
  pjname <- 'all'
  dbpath <- 'all_meta_dir'
  fgpath <- 'all_meta_fig'
  mat_id <- 'NKILC_single'#'all_Imm'#
  # mat_id<-'NK'
  mat_id="NK_nopro"
  mtpath <- 'Mars_batches.txt'
  
  unlink(dbpath, recursive=TRUE)
  if(!dir.exists(dbpath)) dir.create(dbpath)
  scdb_init(dbpath)
  if(!dir.exists(fgpath)) dir.create(fgpath)
  scfigs_init(fgpath)
  
  source("/share/data0/UserData/yangjing/20230710/script/base.R")
  mcell_import_multi_mars(mat_id, mtpath,base_dir="/home/yangjing/SB1921253234/umi_rawdata/",force=TRUE,patch_cell_name=T)
  # mcell_import_multi_mars(mat_id, mtpath,base_dir="/share/data0/UserData/yangjing/SB1921253234/umi_rawdata/NK_SB91",force=TRUE,patch_cell_name=T)
  mat <- scdb_mat(mat_id)
  all=mctoseurat0(mat_id)
  save(all,file="./Seurat/all.rda")
  ##unique cell name
  # colnames(mat@mat)<-add_postfix(colnames(mat@mat))
  # mat@cells<-add_postfix(mat@cells)
  # colnames(mat@ignore_gmat)<-add_postfix(colnames(mat@ignore_gmat))
  # rownames(mat@cell_metadata)<-add_postfix(colnames(mat@ignore_gmat))
  # scdb_add_mat(mat_id, mat)
  # print(dim(mat@mat))
  
  nms = c(rownames(mat@mat), rownames(mat@ignore_gmat)) 
  nms <- unique(c(rownames(mat@mat), rownames(mat@ignore_gmat)))
  # pre_nr_term <- c("^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP")
  # pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
  # pre_ex_genes <- c("MBALAT1", "XIST", "XIST_intron")
  # pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
  # pre_bad_genes
  
  ##remove bad cells  
  genes_RBC <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  cells_RBC <- names(which((apply(mat@mat[intersect(genes_RBC,rownames(mat@mat)),],2,sum))>=1))
  genes_Mt<-foreach(i='^MT-', .combine = c) %do% grep(i, nms, v=T)
  cells_mt= names(which((Matrix::colSums(mat@mat[genes_Mt,])/Matrix::colSums(mat@mat))>0.5))
  mcell_plot_umis_per_cell(mat_id, min_umis_cutoff = 400)
  doublets = read.table("./Seurat/doublet_15.txt",header=T)
  doub_cell<-doublets$cell[which(doublets$doub=="Doublet")]
  cell=c(names(which(Matrix::colSums(mat@mat)>8000)),names(which(Matrix::colSums(mat@mat)<400)))
  unique(c(cells_RBC,cell,cells_mt,doub_cell,cell))->rmcells
  mcell_mat_ignore_cells(mat_id, mat_id, ig_cells = rmcells, reverse = F)
  mat = scdb_mat(mat_id)
  print(dim(mat@mat))
  ##remove bad genes
  nms <- unique(c(rownames(mat@mat), rownames(mat@ignore_gmat)))
  pre_nr_term <- c("^ERCC-","^RPS","^RPL","^MT-","^MTMR","^MTND","^MTRN","^MTCO","^MRPL","^MRPS","^HBA","^HBB","^MTATP")
  pre_nr_genes <- foreach(i=pre_nr_term, .combine = c) %do% grep(i, nms, v=T)
  pre_ex_genes <- c("MALAT1", "XIST", "XIST_intron","Y_RNA")#Y_RNA NKduyou
  pre_bad_genes <- unique(c(pre_nr_genes, pre_ex_genes))
  mcell_mat_ignore_genes(new_mat_id=mat_id, mat_id=mat_id, pre_bad_genes, reverse=F)
  mat = scdb_mat(mat_id)
  print(dim(mat@mat))
  
  genes_anchors = c('FOS','FOSB','NFKBIA','NFKBIZ','JUN','ZFP36','ISG15','HMGB2','STMN1','TOP2A','MKI67','MX1','RSAD2')
    genes_anchors = c('FOS','FOSB','NFKBIA','NFKBIZ','JUN','ZFP36','ISG15','HMGB2','STMN1','TOP2A','MKI67','MX1','RSAD2',
                'CCNB1','CCNB2','CCND2','CCND3','CDKN1A','HMGB1','PCNA')
  tab_fn = "./all_meta_dir/lateral_gmods_2.txt"
  mcell_mat_rpt_cor_anchors(mat_id=mat_id, gene_anchors = genes_anchors, cor_thresh = 0.1,
                            gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = 0.2)
  gcor_mat = read.table('./all_meta_dir/lateral_gmods_2.txt', header=T)
  foc_genes = apply(gcor_mat[, intersect(colnames(gcor_mat),genes_anchors)], 1, which.max)
  
  mcell_add_gene_stat(gstat_id=mat_id, mat_id=mat_id, force=T)
  mcell_gset_filter_varmean(gset_id=mat_id, gstat_id=mat_id, T_vm=0.08, force_new=T)
  mcell_gset_filter_cov(gset_id = mat_id, gstat_id=mat_id, T_tot=100, T_top3=2)
  mcell_plot_gstats(gset_id=mat_id, gstat_id=mat_id)
  gset <- scdb_gset(mat_id)
  
  pst_genes <- names(gset@gene_set)
  pst_nr_term <- c("^AC[0-9]+\\.", "^AP[0-9]+\\.", "^AL[0-9]+\\.", "-AS[0-9]*$", "^MIR[0-9]", "^LINC[0-9]", "^SNOR", "^sno", "^SNR",
                   "^TMSB", "^HIST", "^HSP", "^IFI", "^HLA-", "^ATP", "-", ":", "\\.", '^KIAA',
                   "^IGJ", "^IGH", "^IGK", "^IGL", "^DNAJ", "^GZM", "^CCL", "^XCL", '^FTH', '^FTL', '^LGALS')
  
  pst_nr_genes <- foreach(i=pst_nr_term, .combine = c) %do% grep(i, pst_genes, v=T)
  pst_ex_genes <- c()
  pst_bad_genes <- unique(c(pst_nr_genes, pst_ex_genes, names(foc_genes)))
  pst_add_genes <- c()
  final_genes <- unique(setdiff(pst_genes, pst_bad_genes), pst_add_genes)
  
      # load("/share/data0/UserData/yangjing/20230710/all/all_meta_dir/gset.NKILC_single.Rda")
      # gene=intersect(final_genes,names(object@gene_set))
      # scdb_add_gset(mat_id, gset)
      # gene=intersect(names(gset@gene_set),names(object@gene_set))
      # pst_nr_term_nk<-c("^ZNF","^ZFP","^ZBT","Y_RNA","^TNF","NEAT1","HERC1","Metazoa_SRP","^PRP")
      # pst_nr_genes_nk <- foreach(i=pst_nr_term_nk, .combine = c) %do% grep(i, gene, v=T)
      # final_genes <- setdiff(gene,pst_nr_genes_nk)
  
  gset@gene_set <- rep(1,length(final_genes))
  names(gset@gene_set) <- final_genes
  scdb_add_gset(mat_id,gset)
  
  mcell_add_cgraph_from_mat_bknn(
    mat_id=mat_id,
    gset_id = mat_id,
    graph_id=mat_id,
    K=404,#NK_91=80
    dsamp=T)
  mcell_coclust_from_graph_resamp(
    coc_id=mat_id,
    graph_id=mat_id,
    min_mc_size=100,#50
    p_resamp=0.75, n_resamp=500)
  mcell_mc_from_coclust_balanced(
    coc_id=mat_id,
    mat_id= mat_id,
    mc_id= mat_id,
    K=100, min_mc_size=100, alpha=2)
  
  #### metacel配色
  mc_f <- scdb_mc(mat_id)
  my.col <- c("darkgray","burlywood1","chocolate4","orange","red","purple","blue","darkgoldenrod3","cyan")
  mc_f@colors <- colorRampPalette(my.col)(ncol(mc_f@mc_fp))
  scdb_add_mc(mat_id, mc_f)
  mc_f <- scdb_mc(mat_id)
  
  #### Markers热图(热图显示的log2(lfp)值)
  set_param("scm_mc_mark_k_per_clust",20,"metacell") # default: 5
  mcell_gset_from_mc_markers(gset_id = "pbmc_markers", mc_id = "pbmc")
  mcell_gset_from_mc_markers(gset_id = "NK_markers1", mc_id = mat_id)
  # 单细胞水平
  mcell_mc_plot_marks(mc_id = mat_id, gset_id = "NK_markers1", mat_id = mat_id, plot_cells = F)
  # metacell水平
  mcell_mc_plot_marks(mc_id = "pbmc", gset_id = "pbmc_markers", mat_id = "pbmc", plot_cells = F)
  
  mc<-scdb_mc(mat_id)
  meta <- read.table(mtpath, header=T, row.names = 1, sep="\t")
  meta <- meta[intersect(rownames(meta), rownames(mc@n_bc[rowSums(mc@n_bc) > 0,])),]
  meta <- meta[order(meta$Batch.Set.ID),]
  Total <- rowSums(mc@n_bc[rownames(meta),])
  write.table(cbind(meta,mc@n_bc[rownames(meta),],Total), file.path(fgpath,paste0(pjname,"_nbc.csv")), sep=",", col.names=NA)
  
  
  #### 2D图展示Cells与MCs
  #download.file("http://www.wisdom.weizmann.ac.il/~arnau/metacell_data/Amphimedon_adult/config.yaml","config.yaml")
  #tgconfig::override_params("config.yaml","metacell") #使用这个就不用下面的2-3行代码设置参数
  set_param("mcell_mc2d_K",20,"metacell") # default: 20
  set_param("mcell_mc2d_T_edge",0.017,"metacell") # default: 0.05  all:0.017
  set_param("mcell_mc2d_max_confu_deg",3.5,"metacell") # default: 5 all:3.5
  set_param("mcell_mc2d_edge_asym",FALSE,"metacell") # default: TRUE
  set_param("mcell_mc2d_proj_blur",0.015,"metacell") # default: 0.02
  
  mcell_mc2d_force_knn(mc2d_id=mat_id, mc_id=mat_id, graph_id=mat_id)
  tgconfig::set_param("mcell_mc2d_height", 500, "metacell")
  tgconfig::set_param("mcell_mc2d_width", 500, "metacell")
  mcell_mc2d_plot(mc2d_id = mat_id,sc_cex=0.5)
  
  mc2d <- scdb_mc2d(mat_id)
  
  mc_hc <- mcell_mc_hclust_confu(mc_id=mat_id, graph_id=mat_id)
  mc_sup <- mcell_mc_hierarchy(mc_id=mat_id, mc_hc=mc_hc, T_gap=0.04)
  mcell_mc_plot_hierarchy(mc_id=mat_id, 
                          graph_id=mat_id, 
                          mc_order=mc_hc$order, 
                          sup_mc = mc_sup, 
                          width=3000, heigh=2000, 
                          min_nmc=2, show_mc_ids = T)
  #mcell_mc_export_tab(mc_id = mat_id, gstat_id = "test", mat_id = mat_id, T_fold=2)
  mc<-scdb_mc(mat_id)
  write.table(as.matrix(mc2d@mc_x),"./result_txt/mc_x.txt")
  write.table(as.matrix(mc2d@mc_y),"./result_txt/mc_y.txt")
  write.table(as.matrix(mc2d@sc_x),"./result_txt/sc_x.txt")
  write.table(as.matrix(mc2d@sc_y),"./result_txt/sc_y.txt")
  write.table(mc@mc,"./result_txt/mcmc.txt")
  write.table(mc@annots,"./result_txt/annots.txt")
  write.table(mc@mc_fp,"./result_txt/mc_fp.txt")
  write.table(mc@n_bc,"./result_txt/mcnbc.txt")
  
  lfp = round(log2(mc@mc_fp),2)
  write.table(lfp,"./result_txt/lfp_NK.txt")

#Figure1 
  setwd("./")
  library(Seurat)
  library(sctransform)
  library(ggplot2)
  library(ggsci)
  library(dplyr)
  library(openxlsx)
  library(stringr)
  library(paletteer)
  id="T_single5"#B_single,NK_single,Myloid_single,immu_cell
  source("./script/base.R")
  col_T=paletteer_d("ggsci::category20_d3",13)
  mc2su<-mctoseurat("T_single",id,id)
  save(mc2su,"./mc2su_T.rda")
  #################################metadata load 
  # load("./Seurat/mc2su_T.rda")
  plateinfo=read.table("../plateinfo.txt",sep="\t",header=T)
  # mcinfo=read.table("./Seurat/mcinfo_sure.txt",sep="\t",header=T)
  mcinfo=read.table("./Seurat/mcinfo_T_new.txt",sep="\t",header=T)
  # mc2su$cluster0=mcinfo$annotation[match(mc2su$mc,mcinfo$mc)]
  mc2su$cluster0=mcinfo$cluster0[match(mc2su$mc,mcinfo$mc)];mc2su$cluster1=mcinfo$cluster1[match(mc2su$mc,mcinfo$mc)];
  mc2su@meta.data$sample<-mc2su@meta.data$amp_batch_id
  mc2su@meta.data$cell=colnames(mc2su)
  mc2su@meta.data$batches=plateinfo$seq_batch[match(mc2su$sample,plateinfo$Amp_batch)]
  mc2su@meta.data$Amp_batch=mc2su@meta.data$sample
  mc2su@meta.data$sample=plateinfo$sample.ID[match(mc2su$Amp_batch,plateinfo$Amp_batch)]	
  mc2su@meta.data$id=plateinfo$id[match(mc2su$Amp_batch,plateinfo$Amp_batch)]					
  mc2su@meta.data$day_in_cycle=plateinfo$day.in.cycle[match(mc2su$Amp_batch,plateinfo$Amp_batch)]
  mc2su$disease=plateinfo$type..disease.[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$type_intrauterine=plateinfo$type..intrauterine.[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$ending=plateinfo$ending[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$period_in_cycle=plateinfo$phase[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$Age=plateinfo$age[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$BMI=plateinfo$BMI[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$sorting=plateinfo$sorting.gate[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$stage=plateinfo$stage[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$CD138=plateinfo$CD138.examination[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$step=plateinfo$step[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$FSH=plateinfo$FSH[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$LH=plateinfo$LH[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$E2=plateinfo$E2[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$HCG=plateinfo$HCG[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su$thickness=plateinfo$endometrial.thickness..mm.[match(mc2su$Amp_batch,plateinfo$Amp_batch)] 
  mc2su@meta.data$id=plateinfo$id[match(mc2su$Amp_batch,plateinfo$Amp_batch)]
  mc2su@meta.data$CD45.cell.count=plateinfo$CD45.cell.count[match(mc2su$Amp_batch,plateinfo$Amp_batch)]
  Idents(mc2su)<-mc2su$cluster1;mc2su<-subset(mc2su,idents=c("delete"),invert=T)#doubletfinder doublet
  mc2su$cluster1=factor(mc2su$cluster1,levels=c("naïve T","TEM/Tfh-like","TEM/Th17-like","TM CD4","Treg","ISG15+ CD8","TNFAIP3+ CD8","CD8","CD8p","CD8+ NKT","NKT","cytotoxic T","MAIT"))
  # mc2su$cluster1=factor(mc2su$cluster1,levels=c("NK1","NEAT1+ NK1","SAT1+ NK1","NK2","GNLY+ NK2","NK3","ISG15+ NK","CRTAM+ NK","ILC2","ILC3"))
  # mc2su$cluster1=factor(mc2su$cluster1,levels=c("clasical monocyte","non-clasical monocyte","CD163+ macrophage","DAB2+ macrophage","LYVE1+ macrophage","c3+ macrophage","MARCO+ macrophage","CXCL8 macrophage","osteoclast-like macrophage","cDC1","cDC2","cycling DC","pDC","Granulocytes"))
  # mc2su$cluster1=factor(mc2su$cluster1,levels=c("memory B cell","naive B cell","GC B like","CD83+ B cell","plasma cell"))
  save(mc2su,file="./Seurat/mc2su_T.rda")
  #######################################dimplot of cluster0##########
  load("./Seurat/mc2su_final.rda")
  load("./Seurat/mc2su_myeloid.rda")
  load("./Seurat/mc2su_T.rda")
  load("./Seurat/mc2su_NKILC2.rda")
  load("./Seurat/mc2su_B.rda")
  library(SCpubr)
  library(paletteer)
  Idents(mc2su)<-mc2su$step;mc2su<-subset(mc2su,idents=c("3"),invert=T)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su<-subset(mc2su,idents=c("na"),invert=T)
  Idents(mc2su)<-mc2su$sorting;mc2su<-subset(mc2su,idents=c("CD45+"))
  Idents(mc2su)<-mc2su$disease;mc2su<-subset(mc2su,idents=c("other"),invert=T)

  sc=mc2su@meta.data;sc$x=mc2su@reductions$metacell@cell.embeddings[,1];sc$y=mc2su@reductions$metacell@cell.embeddings[,2]
  p1 = ggplot(sc)+
    geom_point(aes(x = x,y = y,fill = cluster1),shape = 21,color = "grey20",size = 3,stroke = 0.05)+
    theme_bw()+
    xlab("")+ylab("")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank())+
    scale_fill_manual(values = col_B)
  ggsave("T_cluster1_2d1.pdf",p1,width = 13,height = 10,dpi = 300) # save the image

#FigS1
  library(Seurat)
  setwd("./Seurat/")
  load("mc2su_final.rda")
  Idents(mc2su)<-mc2su$step;mc2su<-subset(mc2su,idents=c("3"),invert=T)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su<-subset(mc2su,idents=c("na"),invert=T)
  Idents(mc2su)<-mc2su$sorting;mc2su<-subset(mc2su,idents=c("CD45+"))
  na_samples <- which(is.na(mc2su$cluster1));mc2su <- subset(mc2su, cells = -na_samples)
  Idents(mc2su)<-mc2su$disease;mc2su<-subset(mc2su,idents=c("CE","healthy"))
  pdf("../fig_result/Number_of_UMI_count_per_cell.pdf",height=5,width=5)
  h = hist(mc2su$nCount_RNA,ylab="Cell frequency",xlab="UMI count" ,
      breaks=seq(0, max(mc2su$nCount_RNA) + 100, by=100),col="#69AB32",xlim=c(0,20000),main="Number of UMI count per cell")
  abline(v=median(mc2su$nCount_RNA), col='red', lty=2)
  text(max(mc2su$nCount_RNA), max(h$count)/1, "mean=972", pos=2, col='red')
  dev.off()
  pdf("../fig_result/Number_of_detected_genes_per_cell.pdf",height=5,width=5)
  h = hist(mc2su$nFeature_RNA,ylab="Cell frequency",xlab="Number of detected genes" ,
       breaks=seq(0, max(mc2su$nFeature_RNA) + 50, by=50),col="#F0E356",xlim=c(0,6000),main="Number of detected genes per cell")
  abline(v=median(mc2su$nFeature_RNA), col='red', lty=2)
  text(max(mc2su$nFeature_RNA), max(h$count)/1, "mean=687", pos=2, col='red')
  dev.off()

  library(Seurat)
  setwd("/home/yangjing/20230710/all/Seurat/")
  load("mc2su_final.rda")
  Idents(mc2su)<-mc2su$step;mc2su<-subset(mc2su,idents=c("3"),invert=T)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su<-subset(mc2su,idents=c("na"),invert=T)
  Idents(mc2su)<-mc2su$sorting;mc2su<-subset(mc2su,idents=c("CD45+"))
  Idents(mc2su)<-mc2su$disease;health<-subset(mc2su,idents=c("healthy"))
  mc2samp=as.matrix(table(health$mc,health$sample))
  stage=health$stage[match(colnames(mc2samp),health$sample)]
  stage=gsub("-","_",stage)
  cor1=cor(mc2samp,method="spearman")
  cor_list=list("Early_pro"=c(),"Mid_pro"=c(),"Late_pro"=c(),"Early_sec"=c(),"Mid_sec"=c(),"Late_sec"=c(),"Late_Late_sec"=c(),"menstrual"=c(),"Among_stage"=c())
  for(i in 1:nrow(cor1)){
      tmp=min(i+1,ncol(cor1))
      for(j in tmp:ncol(cor1)){
          if(stage[i]=="Early_pro"&stage[j]=="Early_pro") {
              cor_list$Early_pro=c(cor_list$Early_pro,cor1[i,j])
          }else if (stage[i]=="Mid_pro"&stage[j]=="Mid_pro") {
              cor_list$Mid_pro=c(cor_list$Mid_pro,cor1[i,j])
          }else if (stage[i]=="Late_pro"&stage[j]=="Late_pro") {
              cor_list$Late_pro=c(cor_list$Late_pro,cor1[i,j])}
          else if (stage[i]=="Early_sec"&stage[j]=="Early_sec") {
              cor_list$Early_sec=c(cor_list$Early_sec,cor1[i,j])}
          else if (stage[i]=="Mid_sec"&stage[j]=="Mid_sec") {
              cor_list$Mid_sec=c(cor_list$Mid_sec,cor1[i,j])}
          else if (stage[i]=="Late_sec"&stage[j]=="Late_sec") {
              cor_list$Late_sec=c(cor_list$Late_sec,cor1[i,j])}
          else if (stage[i]=="Late_Late_sec"&stage[j]=="Late_Late_sec") {
              cor_list$Late_Late_sec=c(cor_list$Late_Late_sec,cor1[i,j])}
          else if (stage[i]=="menstrual"&stage[j]=="menstrual") {
              cor_list$menstrual=c(cor_list$menstrual,cor1[i,j])}
          else {
              cor_list$Among_stage=c(cor_list$Among_stage,cor1[i,j])
          }
      }
  }
  cor_list$menstrual=cor_list$menstrual[-4]
  tmp=data.frame("cor"=unlist(cor_list),"type"=rep(names(cor_list),c(21,28,6,21,28,78,3,3,1190)))
  tmp$type=factor(tmp$type,levels=c("Early_pro","Mid_pro","Late_pro","Early_sec","Mid_sec","Late_sec","Late_Late_sec","menstrual","Among_stage"))
  
  library(ggplot2)
  library(ggsignif)
  library(paletteer)
  col_stage=paletteer_d("basetheme::clean",8)
  P1 <- ggplot(tmp,aes(x=type,y=cor,fill=type))+ #”fill=“设置填充颜色
     stat_boxplot(geom = "errorbar",width=0.15)+ #由于自带的箱形图没有胡须末端没有短横线，使用误差条的方式补上
     geom_boxplot(size=0.5,outlier.fill="white",outlier.color="white")+ #size设置箱线图的边框线和胡须的线宽度，fill设置填充颜色，outlier.fill和outlier.color设置异常点的属 +   geom_jitter(aes(fill=type),width =0.2,shape = 21,size=3)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
     #geom_jitter(aes(fill=type),width =0.2,shape = 21,size=3)+ #设置为向水平方向抖动的散点图，width指定了向水平方向抖动，不改变纵轴的值
     scale_fill_manual(values=c(col_stage,"black"))+  #设置填充的颜色
     geom_signif(comparisons = list(c("Among_stage", "Early_pro"),
                                   c("Early_pro", "Mid_pro"),c("Late_pro", "Among_stage")),y_position=c(0.92, 0.87,0.82),map_signif_level=TRUE)+
     theme_classic()+ #背景变为白色
     theme(#legend.position="none", #不需要图例
           axis.text.x=element_text(family="Times",size=12), #设置x轴刻度标签的字体属性
           axis.text.y=element_text(family="Times",size=12,face="plain"), #设置x轴刻度标签的字体属性
           axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴的标题的字体属性
           axis.title.x=element_text(family="Times",size = 15,face="plain"), #设置x轴的标题的字体属性
           plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
           panel.grid.major = element_blank(), #不显示网格线
           panel.grid.minor = element_blank(),
           strip.background=element_rect(color="black"))+
     theme(axis.text.x = element_text(angle = 45, hjust = 1))+
     ylab("Spearman correlation")+xlab("stage")
    ggsave("../fig_result/Spearman correlation.pdf",P1,width = 8, height = 8)

#
  mens_nbc<-read.csv("./mens_cycle_meta_fig/mens_nbc.csv",header = T)
  library(ComplexHeatmap)
  #mens_nbc_4546=mens_nbc[match(unique(umi_5$Amp_batch[which(umi_5$batches=="SB45"|umi_5$batches=="SB46")]),mens_nbc$Batch.Set.ID),];
  #cor_nbc=cor(t(mens_nbc_4546[,6:547]))#45,46 plate#
  colnames(cor_nbc)=mens_nbc$Batch.Set.ID;rownames(cor_nbc)=mens_nbc$Batch.Set.ID;
  plateinfo=read.table("plateinfo.txt",sep="\t",header=T)
  plateinfo<-plateinfo[match(mens_nbc$Batch.Set.ID,plateinfo$Amp_batch),]
  cor_nbc<-cor_nbc[match(plateinfo$Amp_batch[order(plateinfo$id)],rownames(cor_nbc)),match(plateinfo$Amp_batch[order(plateinfo$id)],colnames(cor_nbc))]
  
  tmp1=names(table(plateinfo$id)[which(table(plateinfo$id)>9)])
  plateinfo1=plateinfo[which(is.na(match(plateinfo$id,tmp1))==F),]
  cor_nbc1<-cor_nbc[match(plateinfo1$Amp_batch[order(plateinfo1$id)],rownames(cor_nbc)),match(plateinfo1$Amp_batch[order(plateinfo1$id)],colnames(cor_nbc))]
  
  split = plateinfo$id[match(rownames(cor_nbc),plateinfo$Amp_batch)]
  tmp1=unique(split)
  text = func.vector2List(tmp1)
  names(text) = unique(split)
  row_right=rowAnnotation(Type = anno_textbox(split, text,side="right"))
  pdf(file.path(fgpath,paste0("ziji","_amp_corbetween_sample.pdf")),height=10,width=10)
  pdf(file.path(fgpath,paste0("all","_amp_corbetween_sample.pdf")),height=25,width=25)
  Heatmap(cor_nbc,
    cluster_rows=F,cluster_columns=F,
      show_row_names=F,show_column_names=F,
      width=15,height=15,right_annotation=row_right)
  dev.off()

#Fig2
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(reshape2)
  
  setwd("/share/data0/UserData/yangjing/20230710/all/")
  load("./Seurat/mc2su_NKILC2.rda")

  library(tidyr)
  library(reshape2)

  Idents(mc2su)<-mc2su$step;mc2su<-subset(mc2su,idents=c("3"),invert=T)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su<-subset(mc2su,idents=c("na"),invert=T)
  Idents(mc2su)<-mc2su$sorting;mc2su<-subset(mc2su,idents=c("CD45+"))
  mc2su$stage[which(mc2su$disease%in%c("decb","decp"))]="dec";mc2su$disease[which(mc2su$disease%in%c("decb","decp"))]="healthy"
  Idents(mc2su)<-mc2su$disease;mc2su=subset(mc2su,idents=c("healthy"))
  cell.prop1<-as.data.frame(prop.table(table(mc2su$cluster1,mc2su$sample),2))
  cell.prop2<-as.data.frame(prop.table(table(mc2su$cluster0,mc2su$sample),2))
  cell.prop=rbind(cell.prop1,cell.prop2)
  
  cell.prop=cell.prop[which(cell.prop$Var1 %in% c("NK1","NK2","NK3")),]#"NK Cell",,"ISG15+ NK"
  colnames(cell.prop)<-c("cellTypes","sample","proportion")
  cell.prop$cellTypes=as.character(cell.prop$cellTypes);cell.prop$cellTypes=factor(cell.prop$cellTypes,levels=c("NK1","NK2","NK3"))
  cell.prop$stage=mc2su$stage[match(cell.prop$sample,mc2su$sample)]
  cell.prop$phase=mc2su$period_in_cycle[match(cell.prop$sample,mc2su$sample)]
  cell.prop$stage[which(cell.prop$stage=="Late-Late-sec")]="Premenstrual"
  cell.prop=cell.prop[which(!cell.prop$stage %in% c("menstrual")),]
  cell.prop$stage=factor(cell.prop$stage,levels=c("Early-pro","Mid-pro","Late-pro","Early-sec","Mid-sec","Late-sec","Premenstrual","dec"))
  group_mean <- aggregate(cell.prop$proportion, list(cell.prop$cellTypes,cell.prop$stage), mean)
  
  colnames(group_mean)=c("cellType","stage","proportion")
  group_mean$proportion=group_mean$proportion*100
  nk12ratio=c()
  nk12ratio$stage=levels(group_mean$stage);nk12ratio$cellType=rep("NK1/NK2",8)
  nk12ratio$proportion=as.numeric(group_mean[which(group_mean$cellType %in% c("NK1")),3]/group_mean[which(group_mean$cellType %in% c("NK2")),3])
  nk12=as.data.frame(nk12ratio)
  nk12$stage=factor(nk12$stage,levels=c("Early-pro","Mid-pro","Late-pro","Early-sec","Mid-sec","Late-sec","Premenstrual","dec"))
  
  nk12$proportion[7]=4;nk12$proportion[8]=6.5
  IFNK=read.table("./result_txt/IFNK12ratio.txt",sep="\t",header=T)
  IFNK$menstrual=rep(NA,4)
  rownames(IFNK)=IFNK$X;IFNK=IFNK[,-1];IFNK=as.matrix(t(IFNK));IFNK[,1:3]=IFNK[,1:3]*10
  library(reshape2);group_mean=melt(IFNK);colnames(group_mean)=c("stage","cellType","proportion")
  group_mean$stage=factor(group_mean$stage,levels=c("early.pro","mid.pro","late.pro","early.sec","mid.sec","late.sec","late.late.sec","decidual"))
  col_ratio=c("#80d7e1","#babb72","#bf5046","black")
  col_ratio=c("#9F248FFF","#80d7e1","#babb72","#bf5046","#b781d2","black")
  group_mean[,3]=group_mean[,3]/10;group_mean=rbind(group_mean,nk12)
  p <- ggplot() +
    geom_smooth(data=group_mean,aes(x = stage, y = proportion,color=cellType,group = cellType),size=0.5,method = "loess", se = FALSE)+
    scale_color_manual(values=col_ratio)+
    theme_bw()+
    ggtitle("NK1/NK2 ratio in different stage")
  ggsave("./fig_result/NK1_NK2_ratio1_new1.pdf",p,width=8,height=3)
  ggsave("./fig_result/NK1_NK2_ratio1_IF.pdf",p,width=8,height=3)
  library(patchwork)
  ggsave("./fig_result/NK1_NK2_ratio1_IFcombined.pdf",p/p1,width=8,height=6)

  library(Seurat)
  setwd("/share/data0/UserData/yangjing/20230710/all/")
  load("./Seurat/mc2su_final.rda")
  #load("./Seurat/mc2su_B.rda")
  #load("./Seurat/mc2su_NKILC.rda")
  #load("./Seurat/mc2su_myeloid.rda")
  library(paletteer)
  library(ggalluvial)
  
  col_cluster0=paletteer_d("awtools::spalette",6)
  col_B=paletteer_d("awtools::ppalette",5)
  col_NK<-c("#80d7e1","#babb72","#bf5046","#b781d2","#ece7a3","#f29432","#9c9895","#b38a8f","#6cb25e","#87b2d4")
  col_T=paletteer_d("ggsci::category20_d3",13)
  col_other=paletteer_d("ggsci::springfield_simpsons",14)
  col_Myeloid_cluster0=c("#f9766e","#e1c548","#5fa664","#4e79a6")
  Idents(mc2su)<-mc2su$step;mc2su<-subset(mc2su,idents=c("3"),invert=T)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su<-subset(mc2su,idents=c("na"),invert=T)
  Idents(mc2su)<-mc2su$sorting;mc2su<-subset(mc2su,idents=c("CD45+"))
  Idents(mc2su)<-mc2su$disease;mc2su<-subset(mc2su,idents=c("other"),invert=T)
  mc2su$stage[which(mc2su$stage=="")]=mc2su$disease[which(mc2su$stage=="")]
  Idents(mc2su)<-mc2su$stage;mc2su<-subset(mc2su,idents=c("menstrual","Early-pro","Mid-pro","Late-pro","Early-sec","Mid-sec","Late-sec","Late-Late-sec"))
  
  health=subset(mc2su,idents=c("healthy"))
  health$cluster0=factor(health$cluster0,levels=c("NK Cell","T Cell","monocyte/macrophage/DC","B Cell","ILC","Mast Cell"))
  
  health$stage=ifelse(health$stage=="Late-Late-sec","Premenstrual",health$stage)
  health$stage=factor(health$stage,levels = c("menstrual","Early-pro","Mid-pro","Late-pro","Early-sec","Mid-sec","Late-sec","Premenstrual"))
  
  names(col_cluster0)=levels(health$cluster0)
  names(col_B)=levels(health$cluster1)
  health_all=plot_alluvium(health$stage,health$cluster1,"Menstrual Cycle","Cluster1","Cluster in Healthy",col_B)
  health_all=plot_alluvium(health$stage,health$cluster0,"Menstrual Cycle","Cluster0","Cluster in Healthy",col_cluster0)
  ggsave("./fig_result/Cluster0 in Healthy in B.pdf",health_all,width = 8, height = 8)  
  ggsave("./fig_result/Cluster0 in Healthy.pdf",health_all,width = 8, height = 8)  
  health_all=plot_alluvium(health$stage,health$cluster0,"Menstrual Cycle","Cluster0","Cluster in Healthy",col_Myeloid_cluster0)
  ggsave("/home/yangjing/fig_result/Cluster0 in Healthy in Myeloid_cluster0.pdf",health_all,width = 8, height = 8) 

#Fig3
  library(paletteer)
  library(cowplot)
  library(rstatix)
  library(viridis)
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  setwd("/share/data0/UserData/yangjing/20230710/all/")
  load("./Seurat/mc2su_final.rda")
  load("./Seurat/mc2su_NKILC2.rda")
  load("./Seurat/mc2su_B.rda")
  load("./Seurat/mc2su_T.rda")
  load("./Seurat/mc2su_myeloid.rda")
  
  col_stage=paletteer_d("basetheme::clean",10)
  col_phase=c("#16A84E","#9D33F4","#F5A420")
  col_CD138=c("#6D9FD8", "#E98887")#"#146152","#FF5A33"
  Idents(mc2su)<-mc2su$step;mc2su<-subset(mc2su,idents=c("3"),invert=T)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su<-subset(mc2su,idents=c("na"),invert=T)
  mc2su$stage[which(mc2su$stage=="")]=mc2su$disease[which(mc2su$stage=="")]
  Idents(mc2su)<-mc2su$disease;health=subset(mc2su,idents=c("healthy","CE"))
  df=data.frame()
  sample=unique(mc2su$sample)
  #####单样本相对于另一种疾病状态的所有样本计算roe
  for(i in sample){
    disease=mc2su$disease[match(i,mc2su$sample)]
    cell=colnames(mc2su)[which(mc2su$disease!=disease|mc2su$sample==i)]
    sample_sc=subset(mc2su,cells=cell)
    tbl <- table(sample_sc$cluster0,sample_sc$sample)
    res <- chisq.test(tbl)
    expected = res$expected
    roe <- tbl/expected
    df1 <- roe %>% as.data.frame() %>%
      select(Cluster = Var1, sample = Var2, roe = Freq)
    df_tmp=df1[which(df1$sample==i),]
    # print(df_tmp)
    df=rbind(df,df_tmp)
  }
  
  df$stage=mc2su$stage[match(df$sample,mc2su$sample)]
  df$phase=mc2su$period_in_cycle[match(df$sample,mc2su$sample)]
  df$disease=mc2su$disease[match(df$sample,mc2su$sample)]
  df$disease=factor(df$disease,levels = c("healthy","CE"))
  df1=df[which(is.na(df$roe)==F),]
  a=compare_means(roe~disease,data=df1,group.by ="Cluster")
  P1 <- ggplot(df1,aes(x=disease,y=roe,fill=disease))+ #”fill=“设置填充颜色
      geom_bar(aes(),stat='summary',fun="mean", position=position_dodge(width=0.1),width = 0.7,alpha=0.6)+#,color="black"
      stat_summary(fun.data = "mean_se",geom = "errorbar",colour="black", width = 0.2,position =position_dodge(0.1))+
      # geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),position ="dodge",width=0.2)+
      geom_jitter(aes(color=disease),,width =0.2,shape = 16,size=2)+
      labs(y = "Ro/E", x = "Cluster", color = "disease") +
      geom_hline(aes(yintercept=1), colour="black",linetype="dashed") +
      facet_wrap(~Cluster,ncol = 6,scales="free")+
      scale_fill_manual(values=c(col_CD138))+  #设置填充的颜色
      scale_color_manual(values=c(col_CD138))+ 
      ggtitle("Proportion of cells in different cell type")+ #设置总的标题
      theme_classic()+ #背景变为白色
      stat_compare_means(label="p.signif",label.x=1.5,label.y = 1.8)+
      theme(#legend.position="none", #不需要图例
          axis.text.x=element_text(family="Times",size=12), #设置x轴刻度标签的字体属性
          axis.text.y=element_text(family="Times",size=12,face="plain"), #设置x轴刻度标签的字体属性
          axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴的标题的字体属性
          axis.title.x=element_text(family="Times",size = 15,face="plain"), #设置x轴的标题的字体属性
          plot.title = element_text(family="Times",size=15,face="bold",hjust = 0.5), #设置总标题的字体属性
          panel.grid.major = element_blank(), #不显示网格线
          panel.grid.minor = element_blank(),
          strip.background=element_rect(color="black"))+
      scale_y_continuous(expand = c(0,0))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),panel.spacing = unit(2, "lines"))
    #ylab("Ro/E")+xlab("CD138")
    #ylab("Proportion of cells")+xlab("stage") #设置x轴和y轴的标题
  ggsave("./fig_result/cluster0 diff in CE_Healthy_new.pdf",P1,width = 16, height =5)
  
  library(fmsb)
  library(dplyr)
  library(tidyr)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su=subset(mc2su,idents=c("menstrual phase"),invert=T)
  mc2su$disease_phase=paste0(mc2su$disease,"-",mc2su$period_in_cycle)
  mc2su$disease_phase=factor(mc2su$disease_phase,levels=c("healthy-proliferative phase","healthy-secretory phase","CE-proliferative phase","CE-secretory phase"))
  # mc2su$disease=factor(mc2su$disease,levels=c("healthy","CE"))
  cell.prop<-as.data.frame(prop.table(table(mc2su$cluster0,mc2su$disease_phase),2))
  cell.prop1<-as.data.frame(prop.table(table(mc2su$cluster0,mc2su$disease),2))
  cell.prop=rbind(cell.prop,cell.prop1)
  cell.prop1=cell.prop[which(cell.prop$Var2 %in% c("CE-proliferative phase","CE-secretory phase","CE")),]
  cell.prop2=cell.prop[which(cell.prop$Var2 %in% c("healthy-proliferative phase","healthy-secretory phase","healthy")),]
  cell.prop=c()
  cell.prop$Cluster=cell.prop1$Var1;cell.prop$phase=c(rep("proliferative phase",6),rep("secretory phase",6),rep("all",6));cell.prop$Freq=cell.prop1$Freq/cell.prop2$Freq
  cell.prop=as.data.frame(cell.prop)
  mat <- spread(cell.prop, key =Cluster, value = Freq);rownames(mat)=mat[,1];mat=mat[,-1]
  mat[2,4]=4
  max_values <- rep(4,6)
  min_values <- rep(0,6)
  mat=rbind(max_values,min_values,mat);rownames(mat)[1:2]= c("Max", "Min")
  # mat1=mat*100
  pdf("./fig_result/radarchart_cluster0_disease_phase1.pdf",height=5,width=5)
  radarchart(
    mat, axistype = 1,
    pcol = col_phase[1:3], pfcol = scales::alpha(col_phase[1:3], 0.5), plwd = 2, plty = 1,
    cglcol = "grey", cglty = 1, cglwd = 0.8,pty=32,
    axislabcol = "grey", caxislabels = c(0,1,2,3,4),
    vlcex = 0.7, vlabels = colnames(mat))
  legend(x=0.7, y=1, legend = rownames(mat[-c(1,2),]), bty = "n", pch=20 , col=col_phase[1:3] , text.col = "black", cex=0.5, pt.cex=1.5)
  dev.off()

#Fig5
  setwd("/share/data0/UserData/yangjing/20230710/all/")
  load("./Seurat/mc2su_final.rda")
  library(dplyr)
  library(tidyr)
  library(corrplot)
  Idents(mc2su)<-mc2su$step;mc2su<-subset(mc2su,idents=c("3"),invert=T)
  Idents(mc2su)<-mc2su$period_in_cycle;mc2su<-subset(mc2su,idents=c("na",""),invert=T)
  Idents(mc2su)<-mc2su$sorting;mc2su<-subset(mc2su,idents=c("CD45+"))
  Idents(mc2su)<-mc2su$disease;mc2su=subset(mc2su,idents=c("CE"))
  cell.prop<-as.data.frame(prop.table(table(mc2su$cluster1,mc2su$sample),2))
  mat <- spread(cell.prop, key =Var1, value = Freq)
  rownames(mat)=mat[,1];mat=mat[,-1]
  M <- cor(mat)
  Pval <- cor.mtest(mat)
  M1=M
  M1[which(M1>0.6)]=0.6
  p1=pheatmap::pheatmap(M1,cluster_rows = T,cluster_cols = T,
                    show_rownames = T,show_colnames = T,scale = "none",color = colorRampPalette(colors = c("#80C7F5","white","#E64F40"))(50))#"#FCCCAD"#"#F6CDA8"
  pdf("./fig_result//heatmap_cellprop_correlation_CE2.pdf",height=8,width=8)
  p1
  dev.off()

