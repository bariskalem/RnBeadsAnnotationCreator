########################################################################################################################
## createAnnotationPackage.hg38.R
## created: 2014-02-13
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Annotation package creation for Hg38.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' createAnnotationPackage.hg38
#'
#' Helper function to create annotation package for genome assembly hg38
#' RnBeads annotation for that assembly
#' @return None (invisible \code{NULL}).
#' @author Fabian Mueller
#' @noRd
createAnnotationPackage.hg38 <- function(){

	suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))

	## Genomic sequence and supported chromosomes
	GENOME <- 'BSgenome.Hsapiens.NCBI.GRCh38'
	assign('GENOME', GENOME, .globals)
	CHROMOSOMES <- c(1:22, "X", "Y")
	names(CHROMOSOMES) <- paste0("chr", CHROMOSOMES)
	assign('CHROMOSOMES', CHROMOSOMES, .globals)
	rm(GENOME, CHROMOSOMES)

	## Download SNP annotation
	logger.start("SNP Annotation")
	vcf.files <- paste0(DBSNP.FTP.BASE, "human_9606_b150_GRCh38p7/VCF/00-All.vcf.gz")
	update.annot("snps", "polymorphism information", rnb.update.dbsnp, ftp.files = vcf.files)
	logger.completed()

	## Define genomic regions
	biomart.parameters <- list(
		database.name = "ENSEMBL_MART_ENSEMBL",
		dataset.name = "hsapiens_gene_ensembl",
		db.version = 109,
		required.columns = c(
			"id" = "ensembl_gene_id",
			"chromosome" = "chromosome_name",
			"start" = "start_position",
			"end" = "end_position",
			"strand" = "strand",
			"symbol" = "hgnc_symbol",
			"entrezID" = "entrezgene_id"),
		host = "https://www.ensembl.org")
	logger.start("Region Annotation")
	update.annot("regions", "region annotation", rnb.update.region.annotation,
		biomart.parameters = biomart.parameters)
	rm(biomart.parameters)
	logger.completed()

	## Define genomic sites
	logger.start("Genomic Sites")
	update.annot("sites", "CpG annotation", rnb.update.sites)
	logger.completed()

	## Define MethylationEPIC v2 probe annotations
	logger.start("MethylationEPICv2")
	table.columns <- rnb.get.illumina.annotation.columns("EPICv2")
	update.annot("probesEPICv2", "MethylationEPICv2 annotation", rnb.update.probeEPICv2.annotation,
	             table.columns = table.columns)
	.globals[['sites']][["probesEPICv2"]] <- .globals[['probesEPICv2']][["probes"]]
	logger.completed()
	
	## Add annotation columns to the context probes, showing if they are covered by an assay
	logger.start("Updating Site Annotation with Probes")
	.globals[['sites']] <- rnb.update.site.annotation.with.probes(sites = .globals[['sites']],
	                                                              query.probes = c("probes27", "probes450", "probesEPIC", "probesEPICv2"),
	                                                              platform.names = c("HumanMethylation27", "HumanMethylation450", "MethylationEPIC", "MethylationEPICv2"))
	logger.completed()


	## Create all possible mappings from regions to sites
	logger.start("Mappings")
	update.annot("mappings", "mappings", rnb.create.mappings)
	logger.completed()

	## Export the annotation tables
	rnb.export.annotations.to.data.files()
}
