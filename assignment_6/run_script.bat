@echo off
setlocal enabledelayedexpansion

set BASE_DIR=./results/perturb_search
if not exist "%BASE_DIR%" mkdir "%BASE_DIR%"

REM Fixed parameters
set GENE_EXCHANGE=10
set USE_BRIDGE=0

REM Values to test for perturbation steps
set "PERTURB_LIST=5 10 15"

for %%p in (%PERTURB_LIST%) do (
    set "OUTFILE=%BASE_DIR%/run_perturb_%%p_gene_%GENE_EXCHANGE%_bridge_%USE_BRIDGE%.txt"
    echo Running: n_perturb=%%p, gene_exchange=%GENE_EXCHANGE%, use_bridge=%USE_BRIDGE% --^> !OUTFILE!
    ../assignment_6_grid.exe "!OUTFILE!" %%p %GENE_EXCHANGE% %USE_BRIDGE%
)

echo.
echo All runs completed. Results saved in %BASE_DIR%