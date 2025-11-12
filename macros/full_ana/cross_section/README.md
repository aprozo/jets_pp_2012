# Cross-Section Analysis Macros

This directory contains macros for calculating jet cross-sections by combining matching results and unfolding.

## Files

- **config.h**: Shared configuration file with binning definitions and analysis parameters
  - Contains `pt_reco_bins` and `pt_mc_bins` used by all macros
  - Analysis configuration structure with triggers, jet radii, cuts, etc.

- **cross_section.cpp**: Main analysis macro that combines matching and unfolding
  - Loads unfolding response matrices
  - Applies efficiency corrections (trigger and reconstruction)
  - Calculates differential cross-sections
  - Saves results with covariance matrices

- **plot_cross_section.cpp**: Plotting macro for visualization
  - Plots individual cross-sections
  - Shows unfolding closure tests
  - Compares different triggers

- **makeHistograms.cpp**: Creates 2D histograms (pt vs run) from data

## Usage

### Prerequisites

Make sure you have:
1. **Completed matching analysis** (creates efficiency histograms in `../matching/plots/embedding_root_*/`)
2. **Completed unfolding** (creates response matrices in `../unfolding/response_*.root`)
3. **Conda environment**: `physics` environment with ROOT
4. **RooUnfold**: Installed at `/home/prozorov/install/RooUnfold`

### Environment Setup

The analysis requires both conda ROOT and RooUnfold:

```bash
# Conda provides ROOT
conda activate physics

# RooUnfold provides unfolding capabilities
export ROOUNFOLD_HOME="/home/prozorov/install/RooUnfold"
export ROOT_INCLUDE_PATH="${ROOUNFOLD_HOME}/include:${ROOT_INCLUDE_PATH}"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
```

**Note**: The helper scripts handle environment setup automatically.

### Running the Analysis

#### Option 1: Using helper scripts (recommended)

```bash
# Test single configuration first
./test_single_config.sh JP2 0.2

# If successful, run full analysis (all triggers + jet radii)
./run_cross_section.sh

# Create plots
./run_plots.sh
```

#### Option 2: Manual execution

```bash
# Set up environment
conda activate physics
export ROOUNFOLD_HOME="/home/prozorov/install/RooUnfold"
export ROOT_INCLUDE_PATH="${ROOUNFOLD_HOME}/include:${ROOT_INCLUDE_PATH}"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}"

# Run in ROOT
root -l -b -q \
  -e 'gSystem->Load("libRooUnfold"); gSystem->AddIncludePath(Form("-I%s/include", getenv("ROOUNFOLD_HOME")));' \
  cross_section.cpp
```

Or for a specific configuration:
```bash
root -l
gSystem->Load("libRooUnfold");
.L cross_section.cpp
calculateCrossSection("JP2", "0.6")  // Single trigger and jet R
```

### Creating Plots

```bash
# Using helper script (recommended)
./run_plots.sh

# Or manually
conda activate physics
export ROOUNFOLD_HOME="/home/prozorov/install/RooUnfold"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
root -l -b -q \
  -e 'gSystem->Load("libRooUnfold");' \
  plot_cross_section.cpp
```

## Workflow

The complete analysis workflow is:

1. **Matching Analysis** (`../matching/matching.cpp`)
   - Creates matched MC trees
   - Calculates trigger and reconstruction efficiencies
   - Outputs: `merged_matching_*.root` and efficiency histograms

2. **Unfolding** (`../unfolding/unfold.cxx`)
   - Creates response matrices from matched jets
   - Performs closure tests
   - Outputs: `response_*.root`

3. **Cross-Section Calculation** (`cross_section.cpp`)
   - Unfolds measured spectra
   - Applies efficiency corrections
   - Normalizes to cross-section
   - Outputs: `cross_section_*.root`

4. **Plotting** (`plot_cross_section.cpp`)
   - Visualizes results
   - Outputs: PDF and PNG plots

## Output Files

### Cross-Section Results (`cross_section_<trigger>_R<jetR>.root`)

Contains:
- `unfolded_spectrum`: Unfolded jet pT spectrum
- `cross_section`: Corrected differential cross-section
- `measured_spectrum`: Raw measured spectrum
- `truth_spectrum`: MC truth spectrum
- `trigger_efficiency`: Trigger efficiency vs pT
- `reconstruction_efficiency`: Reconstruction efficiency vs pT
- `covariance_matrix`: Unfolding covariance matrix

## Configuration

Edit `config.h` to modify:
- Binning (`pt_reco_bins`, `pt_mc_bins`)
- Analysis parameters (vertex cuts, iterations, etc.)
- File paths
- Bad run lists

## Notes

- All macros use the shared binning from `config.h`
- The config can be included in other macros: `#include "config.h"`
- Efficiencies are read from matching results
- Response matrices are read from unfolding results
- Results are saved per trigger and jet radius
- Binning: Uses unfolding bins (48 reco + 14 MC) from config.h
- Efficiency interpolation: Efficiency histograms are interpolated to match unfolded spectrum bins

## Troubleshooting

### RooUnfold Library Not Found

**Error**: `Error in <TUnixSystem::DynamicPathName>: libRooUnfold`

**Solution**: Ensure RooUnfold environment variables are set:
```bash
export ROOUNFOLD_HOME="/home/prozorov/install/RooUnfold"
export LD_LIBRARY_PATH="${ROOUNFOLD_HOME}/lib:${ROOTSYS}/lib:${LD_LIBRARY_PATH}"
```

### Missing Input Files

**Error**: `Error: Cannot open response file` or `Cannot open efficiency file`

**Solution**:
1. Check that matching analysis has been completed: `ls ../matching/plots/embedding_root_*/`
2. Check that unfolding has been completed: `ls ../unfolding/response_*.root`
3. Run the prerequisites first:
   ```bash
   cd ../matching && root -l matching.cpp
   cd ../unfolding && root -l unfold.cxx
   ```

### Conda Environment Issues

**Error**: `physics environment not found`

**Solution**: Create the environment:
```bash
conda create -n physics root -c conda-forge
conda activate physics
```

### Cross-Section Values Seem Wrong

**Check**:
1. Verify efficiency corrections are reasonable (0 < eff < 1)
2. Check unfolding closure in plots (unfolded/truth should be ~1)
3. Compare with reference cross-sections if available
4. Check eta acceptance normalization (should be 2*R)

### Test Script Fails

**Solution**: Run test script with single configuration to isolate the issue:
```bash
./test_single_config.sh JP2 0.2
```

Check the output for specific error messages about missing files or failed calculations.
