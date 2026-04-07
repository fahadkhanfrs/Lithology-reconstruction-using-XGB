Applications of rhob as a target:

    Lithostratigraphy & Correlation: It helps identify different rock bodies based on their density characteristics (e.g., distinguishing high-density carbonate or limestone from lower-density shale or porous sandstone).

    Sequence Stratigraphy: RHOB is used alongside Gamma Ray (GR) to define sequence boundaries (SB) and maximum flooding surfaces (MFS) by recognizing changes in sediment coarsening/fining upward sequences.

    Cyclostratigraphy: RHOB is essential for identifying Milankovitch cycles (astronomical cycles) and calculating Sediment Accumulation Rates (SAR) by identifying cyclicity in sedimentation patterns.

Certain logs such as density-derived porosity introduce information leakage when predicting RHOB, so they were excluded to ensure fair evaluation.

Adding additional logs does not always improve performance due to redundancy and noise, especially under depth-based generalization. Simpler feature sets can sometimes perform better.

“Let me sample every 10 meters”

That gave:

Accuracy ≈ 0.77
Step 2 (new idea)

Now you say:

“Instead of sampling blindly, let me ask the model:
WHERE are you confused?”

Then:

Model says: “I’m unsure at depth = 42”
You sample that point
Model improves
Final idea:

“We don’t need many samples, just the right samples”

oh, wow.

Uncertainty-based sampling without geological context does not outperform uniform sampling.

Adaptive sampling requires contextual geological information to be effective; depth-only models fail to capture facies transitions.

Siltstone acts as a transitional facies between sand and mud, and is often misclassified due to overlapping characteristics.