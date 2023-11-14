########################################################################################################################
## probesEPIC.R
## created: 2015-11-06
## creator: Yassen Assenov
## ---------------------------------------------------------------------------------------------------------------------
## Initializes the Infinium MethylationEPIC probe definition tables by loading them from the Illumina website.
########################################################################################################################

## F U N C T I O N S ###################################################################################################

#' rnb.update.probeEPICv2.annotation
#'
#' Creates probe annotation tables for MethylationEPICv2.
#'
#' @param table.columns Expected columns in the probe annotation table, given as a named \code{character} vector.
#' @return \code{list} of two items:
#'         \describe{
#'           \item{\code{"probes"}}{\code{GRangesList} instance containing probe annotations, one \code{GRanges} per
#'                chromosome.}
#'           \item{\code{"controls"}}{\code{data.frame} with control probe annotation.}
#'         }
#' @author Yassen Assenov
#' @noRd
rnb.update.probeEPICv2.annotation <- function(table.columns) {
  
  ## Load manifest file
  man.file <- "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/MethylationEPIC%20v2%20Files.zip"
  
  ## Create temporary directory
  ## TODO: Remove the files after this section is done.
  dir.epic.v2 <- file.path(.globals[['DIR.PACKAGE']], "temp", "EPIC_v2")
  if (file.exists(dir.epic.v2)){
    if (unlink(dir.epic.v2, TRUE, TRUE) != 0) {
      logger.error(c("Could not remove existing temporary directory", dir.epic.v2))
    }
  }
  if (!dir.create(dir.epic.v2, FALSE, recursive = TRUE)){
    logger.error(c("Could not create temporary directory", dir.epic.v2))
  }

  ## Download the manifest file from Illumina website
  destfile.zip <- file.path(dir.epic.v2, "MethylationEPIC%20v2%20Files.zip")
  if (file.exists(destfile.zip)) {
    logger.status(c("File", destfile.zip, "already downloaded"))
  } else {
    if (download.file(man.file, destfile.zip, quiet = TRUE, mode = "wb") != 0) {
      logger.error(c("Could not download", man.file))
    }
    logger.status(c("Downloaded", man.file))
  }
  
  ## Assign destfile to result
  result <- unzip(zipfile=destfile.zip, exdir=dir.epic.v2)[2]
  
	## Identify methylation and control probe annotation tables
	txt <- scan(result, what = character(), sep = "\n", quiet = TRUE)
	assay.start <- grep("^\\[Assay\\]", txt)
	assay.start <- ifelse(length(assay.start) != 0, assay.start[1], 0L)
	controls.start <- grep("^\\[Controls\\]", txt)
	if (!(length(controls.start) == 1 && controls.start > assay.start + 1)) {
		logger.error("Missing or invalid [Controls] section")
	}
	rm(txt); invisible(gc())

	## Load the methylation and control probe annotation tables
	probe.infos <- read.csv(result, skip = assay.start, nrows = controls.start - assay.start - 2, check.names = FALSE,
		stringsAsFactors = FALSE)
	probe.infos <- probe.infos[, sapply(probe.infos, function(x) { !all(is.na(x)) })]
	if (!identical(colnames(probe.infos), names(table.columns))) {
		logger.error("Unexpected columns in the probe definition table")
	}
	colnames(probe.infos) <- table.columns[colnames(probe.infos)]
	control.probe.infos <- read.csv(result, header = FALSE, skip = controls.start, stringsAsFactors = FALSE)
	control.probe.infos <- control.probe.infos[, sapply(control.probe.infos, function(x) { !all(is.na(x)) })]
	logger.status("Loaded the probe definition tables from Illumina's web site")

	## Validate probe.infos columns and some of the values
	probe.infos[, "Design"] <- as.factor(probe.infos[, "Design"])
	probe.infos[, "Color"] <- as.factor(probe.infos[, "Color"])
	probe.infos[, "Next Base"] <- as.factor(probe.infos[, "Next Base"])
	probe.infos[, "CGI Relation"] <- as.factor(probe.infos[, "CGI Relation"])
	probe.infos[, "DMR"] <- as.factor(probe.infos[, "DMR"])
	probe.infos[, "Enhancer"] <- as.factor(probe.infos[, "Enhancer"])
	probe.infos[, "Regulatory Feature Group"] <- as.factor(probe.infos[, "Regulatory Feature Group"])
	probe.infos <- rnb.probes.fix.infinium.columns(probe.infos, platform = "EPICv2")

	## Validate the control probes
	if (ncol(control.probe.infos) != 5) {
		logger.error("Unexpected number of columns in the control probe definition table")
	}
	colnames(control.probe.infos) <- INF.CONTROL.PROBE.TABLE.COLUMNS
	if (anyDuplicated(control.probe.infos$ID) != 0) {
		logger.error("Duplicated IDs in the control probe definition table")
	}
	if (!identical(sort(unique(control.probe.infos[, 2])), unname(RnBeads:::EPIC.CONTROL.TARGETS))) {
		logger.error("Unexpected values for Target in the control probe definition table")
	}
	control.probe.infos[, 2] <- factor(control.probe.infos[, 2], levels = unname(RnBeads:::EPIC.CONTROL.TARGETS))
	control.probe.infos[, 3] <- factor(RnBeads:::capitalize(tolower(control.probe.infos[, 3])))
	control.probe.infos[, 5] <- control.probe.infos[, 5] == "AVG"
	control.probe.infos <- rnb.update.controlsEPICv2.enrich(control.probe.infos)
	logger.status("Processed control probe annotation table")

	## Add information about CpG counts and GC content in the neighborhood, context, overlaps with SNPs
	## TODO: Is "Guessed strand" annotation necessary? It's not being utilized in RnBeads 
	# saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-1.RDS"))
	# probe.infos <- rnb.update.probe.annotation.guess.strand(probe.infos) # GUESSED STRAND
	# i <- which(probe.infos[, "Strand"] != probe.infos[, "Guessed Strand"])
	# if (length(i) != 0) { # GUESSED STRAND
	# 	i <- table(substr(rownames(probe.infos)[i], 1, 2))
	# 	i <- paste0(names(i), " (", i, " probes)")
	# 	logger.warning(paste("The Strand info for the following probe types might be wrong:", i))
	# }
	# rm(i)
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-2.RDS"))
	# probe.infos <- rnb.update.probesEPICv2.snps(probe.infos) # TODO: Is this needed? The information is already present on the manifest
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-3.RDS"))
	probe.infos <- rnb.update.probe.annotation.cpg.context(probe.infos)
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-4.RDS"))
	probe.infos <- rnb.update.probe.annotation.snps(probe.infos)
#	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPIC-5.RDS"))

	## Add data on cross-hybridization
	## No data available for EPIC v2
	# probe.infos[, "Cross-reactive"] <- rnb.update.probe.annotation.cr(probe.infos[, "ID"], "MethylationEPICv2")

	## Remove NA chromosomes
	## TODO: These are chr0 and chrM. What do we do with these?
	probe.infos <- probe.infos[(which(!is.na(probe.infos[, "Chromosome"]))), ]

	## Convert to GRangesList
	saveRDS(probe.infos, file.path(.globals[['DIR.PACKAGE']], "temp", "probesEPICv2-7.RDS"))
	probes.gr <- rnb.probe.infos.to.GRanges(probe.infos)

	return(list(probes = probes.gr, controls = control.probe.infos))
}

########################################################################################################################

#' rnb.update.controlsEPICv2.enrich
#'
#' Extends the given table of control probe annotations by adding (or replacing) four columns and filling them with
#' values depending on the values of some of the original columns.
#'
#' @param control.probe.infos Table of control probe annotation in the form of a \code{data.frame} containing at least
#'                            the following columns: \code{"Target"}, \code{"Color"} and \code{"Description"}.
#' @return Enriched probe annotation; the given \code{data.frame} with four added or replaced columns:
#'         \code{"Evaluate Green"}, \code{"Evaluate Red"}, \code{"Expected Intensity"} and \code{"Sample-dependent"}.
#' @author Pavlo Lutsik
#' @noRd
rnb.update.controlsEPICv2.enrich <- function(control.probe.infos) {

	## Control probe colors associated with the evaluation of the Green channel
	CONTROL.COLORS.GREEN <- c("Black", "Blue", "Cyan", "Green", "Lime", "Limegreen", "Skyblue")

	## Control probe colors associated with the evaluation of the Red channel
	CONTROL.COLORS.RED <- c("Gold", "Orange", "Purple", "Red", "Tomato", "Yellow")

	## Add columns Evaluate Green and Evaluate Red
	control.probe.infos[["Evaluate Green"]] <- "-"
	control.probe.infos[["Evaluate Red"]] <- "-"
	i <- grep("^DNP", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("^Biotin", control.probe.infos[, "Description"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Color"]] %in% CONTROL.COLORS.GREEN)
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "-"
	i <- (control.probe.infos[["Color"]] %in% CONTROL.COLORS.RED)
	control.probe.infos[i, "Evaluate Green"] <- "-"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	i <- grep("^NEGATIVE", control.probe.infos[, "Target"])
	control.probe.infos[i, "Evaluate Green"] <- "+"
	control.probe.infos[i, "Evaluate Red"] <- "+"
	control.probe.infos[["Evaluate Green"]] <- factor(control.probe.infos[["Evaluate Green"]], levels = c("-", "+"))
	control.probe.infos[["Evaluate Red"]] <- factor(control.probe.infos[["Evaluate Red"]], levels = c("-", "+"))

	## Add column Expected Intensity
	control.probe.infos[["Expected Intensity"]] <- as.character(NA)
	i <- control.probe.infos[, "Target"] %in% c("NEGATIVE", "TARGET REMOVAL", "RESTORATION")
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- c("BISULFITE CONVERSION II", "SPECIFICITY II", "EXTENSION", "NON-POLYMORPHIC")
	i <- control.probe.infos[, "Target"] %in% i
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- control.probe.infos[, "Target"] %in% paste("NORM", c("A", "C", "G", "T"), sep = "_")
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\((High)|(20K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("\\((Medium)|(5K)\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Medium"
	i <- grep("\\(Low\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Low"

	#i <- grep("\\((Bkg)|(5K)\\)$", control.probe.infos[, "Description"])
	i <- grep("\\(Bkg\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^BS Conversion I[- ]C", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^BS Conversion I[- ]U", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	i <- grep("^GT Mismatch.+\\(PM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "High"
	i <- grep("^GT Mismatch.+\\(MM\\)$", control.probe.infos[, "Description"])
	control.probe.infos[i, "Expected Intensity"] <- "Background"
	control.probe.infos[["Expected Intensity"]] <- factor(control.probe.infos[["Expected Intensity"]])

	## Add column Sample-dependent
	control.probe.infos[["Sample-dependent"]] <-
		!(control.probe.infos[["Target"]] %in% RnBeads:::CONTROL.TARGETS.SAMPLE.INDEPENDENT)

	## Add column Index
	control.probe.infos[["Index"]][order(control.probe.infos$Target)] <-
		unlist(sapply(sort(unique(control.probe.infos$Target)), function(target) {
				1:length(which(control.probe.infos$Target==target))
			}
		))

	return(control.probe.infos)
}

########################################################################################################################

#' rnb.update.probe.annotation.guess.strand
#' 
#' Updates the MethylationEPIC probe annotation table by adding a column named \code{Guessed Strand} and setting it to
#' the best guess based on the probe sequence and design type.
#' 
#' @param probe.infos Probe annotation table for MethylationEPIC v2 in the form of a \code{data.frame}.
#' @return The updated probe annotation table.
#' @author Yassen Assenov
#' @noRd
rnb.update.probe.annotation.guess.strand <- function(probe.infos) {
	genome.data <- rnb.genome.data() # Gets the targeted genome assembly sequence.
	# print(genome.data)
	guessed <- data.frame(
		"Guessed Strand" = factor(rep("*", nrow(probe.infos)), levels = c("+", "-", "*")),
		"Mismatches A" = 0L, "Mismatches B" = 0L, check.names = FALSE)

	for (chromosome in names(.globals[['CHROMOSOMES']])) {
		chrom.sequence <- genome.data[[chromosome]]
		for (pr.design in c("I", "II")) {
			i <- which(probe.infos[["Chromosome"]] == chromosome & probe.infos[["Design"]] == pr.design)
			if (length(i) != 0) {
				alleles.A <- as.character(probe.infos[i, "AlleleA Probe Sequence"])
				if (pr.design == "I") {
					alleles.B <- as.character(probe.infos[i, "AlleleB Probe Sequence"])
				} else {
					alleles.B <- NULL
				}
				loci <- probe.infos[i, "Location"]
				guessed[i, ] <- rnb.seq.guess.strands(chrom.sequence, loci, pr.design, alleles.A, alleles.B)
				rm(alleles.A, alleles.B, loci)
			}
		}
	}
	rm(chromosome, chrom.sequence, pr.design, i)

	i <- which(probe.infos$Context == "Other")
	guessed[i, "Guessed Strand"] <- "*"
	guessed[i, "Mismatches A"] <- NA
	guessed[i, "Mismatches B"] <- NA
	cbind(probe.infos, guessed)
}

########################################################################################################################

#' rnb.update.probesEPICv2.snps
#' 
#' Sets chromosome and location information for the SNP probes in MethylationEPIC v2 array, copying it from Infinium EPIC v1.
#' 
#' @param probe.infos Probe annotation table for MethylationEPIC v2 in the form of a \code{data.frame}.
#' @return The updated probe annotation table.
#' @author Yassen Assenov
#' @noRd
rnb.update.probesEPICv2.snps <- function(probe.infos) {
	## TODO: For EPIC V2
	## Identify the SNP probes in EPIC v1 = EPICv1_loci
	probesEPIC.rs <- rnb.annotation2data.frame(.globals$sites$probesEPIC)
	probesEPIC.rs <- probesEPIC.rs[grep("^rs", rownames(probesEPIC.rs)), ]
	extended.ids.epic <- paste(probesEPIC.rs$ID, probesEPIC.rs$Design, probesEPIC.rs$Color, sep = ".")

	## Use the MAF data from the manifest 
	probe.infos[["SNPs MAF"]] <- sapply(seq_len(nrow(probe.infos)), function(i) {
		SNP_DISTANCE <- as.numeric(strsplit(as.character(probe.infos[i, "SNP Distance"]), ";")[[1]])
		SNP_MAF <- as.numeric(strsplit(as.character(probe.infos[i, "SNP Minor Allele Frequency"]), ";")[[1]])
		
		sum(SNP_DISTANCE < 3 & SNP_MAF > 0.05)
	})
	
	## Map the SNP probes in EPIC v2 to the corresponding SNP probes in EPIC v2
	i <- grep("^rs", probe.infos$ID)
	extended.ids.epicv2 <- paste(probe.infos$ID[i], probe.infos$Design[i], probe.infos$Color[i], sep = ".")
	if (!all(extended.ids.epicv2 %in% extended.ids.epic)) {
		logger.error("Not all SNP probes in EPIC v2 found in EPIC v1")
	}
	i850k <- unname(sapply(extended.ids.epicv2, function(x) { which(extended.ids.epic == x) }))
	
	## Transfer the information about the SNP probes to EPIC v2
	probe.infos[i, "Genome Build"] <- 38L
	probe.infos[i, "Chromosome"] <- probesEPIC.rs[i850k, "Chromosome"]
	probe.infos[i, "Location"] <- probesEPIC.rs[i850k, "Start"]
	probe.infos[i, "Strand"] <- probesEPIC.rs[i850k, "Strand"]
	probe.infos[i, "MethylationEPIC"] <- TRUE
	probe.infos
}
