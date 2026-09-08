# Experimental context: dissected E14.5 mouse MGE

The PCDH19 workflow analyzed in the HiCAT pilot uses **dissected mouse medial ganglionic eminence (MGE), embryonic day E14.5**. These samples are **not organoids**. The user explicitly confirmed this context during the Step 07 hierarchy-validation review on 2026-09-08. Repository, directory, and environment names containing `organoid` are legacy infrastructure names and do not describe this experiment.

The expanded pilot contains 12,000 cells: 1,000 from each of 12 registered samples, retaining 19,071 measured genes. Its counts originate from the approved Step 02 QC-filtered dataset. Step 06 supplies existing unintegrated coordinates for display. Step 07 does not infer culture age, adult interneuron subtype, or anatomical location from RNA alone.

The user supplied a textual provisional annotation reference. It is a set of hypotheses for independent review, not labels for fitting clusters, tuning HiCAT parameters, or defining ground truth. Developmental identity must be distinguished from cell cycle, maturation, stress, sample association, and non-neural cell populations. Adult marker absence at E14.5 is not sufficient to reject an immature lineage. In particular, basal/oRG-like gene expression does not establish anatomical oRG identity in embryonic mouse MGE.

Step 07 uses the existing pilot, existing second-seed result, one controlled fine-threshold comparison within the fixed baseline coarse parents, and full-data marker/display checks. It does not cluster the full dataset, transfer pilot labels to full-data cells, lock annotations, remove populations, integrate batches, or regress biological programs.
