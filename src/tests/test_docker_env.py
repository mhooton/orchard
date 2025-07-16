import os
import subprocess
import sys

print("=== Docker Environment Test ===")
print(f"Python executable: {sys.executable}")
print(f"ORCHARD_PATH: {os.environ.get('ORCHARD_PATH', 'NOT SET')}")
print(f"Working directory: {os.getcwd()}")

# Test casutools
try:
    result = subprocess.run(['which', 'wcsfit'], capture_output=True, text=True)
    print(f"wcsfit location: {result.stdout.strip()}")
except Exception as e:
    print(f"wcsfit test failed: {e}")

# Test data access
data_path = "/data/SPECULOOSPipeline"
if os.path.exists(data_path):
    print(f"Data directory accessible: {os.listdir(data_path)}")
else:
    print("Data directory NOT accessible")

print("=== Test Complete ===")