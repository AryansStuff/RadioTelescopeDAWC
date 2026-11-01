#!/usr/bin/env bash
set -e
sudo apt update
sudo xargs -a apt-packages.txt apt install -y

echo 'blacklist dvd_usb_rtl28xxu' | sudo tee /etc/modprobe.d/blacklist-rtl.conf

python3 -m venv --system-site-packages .env
.env/bin/pip install --upgrade pip setuptools
.env/bin/pip install -r requirements.txt

echo "Done, run using .env/bin/python recordspectrum.py, then use runthis.py"
