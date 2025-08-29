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


