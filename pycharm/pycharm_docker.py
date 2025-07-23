import sys
import os

print(f"Working directory: {os.getcwd()}")
print(f"Environment PYTHONPATH: {os.environ.get('PYTHONPATH', 'Not set')}")
print(f"sys.path: {[p for p in sys.path if 'project' in p]}")

try:
    from photometry import casutools
    print("✅ Pipeline modules accessible!")
except ImportError as e:
    print(f"❌ Import failed: {e}")
    print("Available directories in /opt/project:")
    print(os.listdir('/opt/project'))
    print("Available directories in /opt/project/src:")
    print(os.listdir('/opt/project/src'))
