import yaml
import sys

config_file = sys.argv[1]
with open(config_file, 'r') as f:
    config = yaml.safe_load(f)

for k, v in config.items():
    if isinstance(v, dict):
        for sub_k, sub_v in v.items():
            if isinstance(sub_v, list):
                sub_v_str = " ".join(map(str, sub_v))
                print(f"{k}_{sub_k}='{sub_v_str}'")
            else:
                print(f"{k}_{sub_k}='{sub_v}'")
    else:
        print(f"{k}='{v}'")