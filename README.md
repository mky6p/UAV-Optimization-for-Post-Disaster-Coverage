# 1. Process JSON data
matlab -r "run('BuildingProcessor.m')"

# 2. Generate initial UAV positions (e.g., 5 UAVs)
matlab -r "run('Signal_Strength.m')"  # Enter 5 when prompted

# 3. Optimize UAV positions
matlab -r "run('Algorithm.m')"

# 4. Run sensitivity analysis
matlab -r "run('PositionSensitivityAnalysis.m')"
