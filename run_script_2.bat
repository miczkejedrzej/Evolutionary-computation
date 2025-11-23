@echo off
setlocal enabledelayedexpansion

set BASE_DIR=assignment_6\results\perturb_search
if not exist "%BASE_DIR%" mkdir "%BASE_DIR%"

REM Fixed parameters
set PERTURB=10
set USE_BRIDGE=0

REM Values to test for perturbation steps
set "GENE_LIST=5 10 15"

for %%g in (%GENE_LIST%) do (
    set "OUTFILE=%BASE_DIR%/run_perturb_%PERTURB%_gene_%%g_bridge_%USE_BRIDGE%.txt"
    echo Running: n_perturb=%PERTURB%, gene_exchange=%%g, use_bridge=%USE_BRIDGE% --^> !OUTFILE!
    .\assignment_6_grid.exe "!OUTFILE!" %PERTURB% %%g %USE_BRIDGE%
)

echo.
echo All runs completed. Results saved in %BASE_DIR%