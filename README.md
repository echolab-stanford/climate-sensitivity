# climate-sensitivity

Replication repository for the study: *"No consistent evidence society is becoming less sensitive to climate"*

---

## 1. System Requirements

### Software Dependencies
- **Tested on:**
     - R version 4.3.1
- **R packages** (automatically installed via the scripts if not present)

### Operating System Compatibility
- **Tested on:**
  - macOS Sequoia 15.3.2

### Hardware Requirements
- No special hardware is required. Scripts have been tested on typical desktop configurations:
  - 8+ GB RAM
  - Quad-core CPU
  - 2+ GB free disk space
---

## 2. Installation Guide

### Instructions

1. **Clone the repository**:
   ```
   git clone https://github.com/yourusername/climate-sensitivity.git
   cd climate-sensitivity
   ```
2. **Download the data from Dropbox**:
   [Dropbox Link](https://www.dropbox.com/scl/fo/fi0kets79nq0r7ufai23v/ACsWw4K1R-tvR6oGN5Phs8U?rlkey=j9ft96t315w6xkx0euhpf6fdm&st=njofsh3d&dl=0)
3. **Place data files** inside the `data/` directory:
   ```
      climate-sensitivity/
   ├── data/
   │   ├── [Dropbox files here]
   └── scripts/
   ```
#### Typical install time
- Installation and setup take approxiamtely 5-10 mins on a typical desktop computer internet access.

## 3. Demo
### Instruction to Run Demo

1. After setup is complete, run the following scripts from `scripts/`: 
   - `main_text_figure_1.R`
   - `main_text_figure_2.R`
   - `main_text_figure_3_part_1-4.R`
   - `main_text_figure_4.R`
   - `main_text_figure_5.R`

2. **Outputs** from each script will be saved under:
   - `fig/main/`

### Expected Run Time
- Total runtime: ~1min-5min
  - Individual script run times vary from 2 seconds to 10 seconds

**Note**: The current setup produces **only the main figures**. Supplementary figures and tables will be included in a future update or made available via a separate script batch.

## 4. Instructions for Use
### How to run the software on your own data
1. Replace the contents of the `data/` folder with your dataset following the same file naming conventions used in the Dropbox files
2. Adjust file paths in the R scripts if necessary
3. Then re-run the scripts in the `scripts/` folder as shown in the Demo
   
### How to Replicate Results

To fully reproduce the results from the paper:

1. **Clone the repo** and **download the Dropbox data** as described above.
2. **Ensure your R environment** matches the versions listed under "System Requirements".
3. **Execute the scripts** in the order listed in the Demo section
4. **Compare output figures** to those in the paper's main materials

5. **Download data from Dropbox:**  
   [Dropbox Link](https://www.dropbox.com/scl/fo/fi0kets79nq0r7ufai23v/ACsWw4K1R-tvR6oGN5Phs8U?rlkey=j9ft96t315w6xkx0euhpf6fdm&st=njofsh3d&dl=0)

6. **Place the downloaded Dropbox files** into the same folder as the cloned GitHub repository (i.e., inside `climate-sensitivity/data/`).

7. **Set your working directory** to the repository root (e.g., `~/climate-sensitivity`).

8. **Run the following R scripts in order (from within the `scripts/` folder):**
   - `main_text_figure_1.R`
   - `main_text_figure_2.R`
   - `main_text_figure_3_part_1-4.R`
   - `main_text_figure_4.R`
   - `main_text_figure_5.R`

9. **Outputs** from each script will be saved under:
   - `fig/main/`
 

