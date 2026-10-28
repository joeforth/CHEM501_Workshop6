# CHEM501 - Workshop 6
## Git, GitHub, and FAIR Data
In this workshop we'll explore using Git and GitHub to share code and document your projects. We'll explore how to get our Nicla to send data over WiFi and you'll start planning your projects. 

## Setup
### Required Software
Install [Git](https://git-scm.com/install) - you can find it [here](https://git-scm.com/install). 

If you're using Windows, there's a nice installer that you just double-click.

If you're using a Mac, the easiest way of doing this is running the following command in a Terminal:
```
xcode-select --install
```
This will take 15 minutes or so, so it's recommended you do this before the Workshop or (at the very least) during the break.

### Required Arduino Modules
As well as the Arduino boards and modules you installed for [Workshop 1]((https://github.com/joeforth/CHEM501_Workshop1)), you'll also need to install the following modules / libraries on your Arduino IDE:

1. WiFiNINA
1. ArduinoMqttClient

### Required Python Modules
Only a single Python module (`paho-mqtt`) is required to run the code as is provided to you. Run the following terminal commands to create a new conda environment and install it:

```
conda create -n chem501_w6
conda activate chem501_w6
conda install python paho-mqtt
```

Alternatively, you can install using pip using either:

```
pip install paho-mqtt
```

or:

```
git clone https://github.com/joeforth/CHEM501_Workshop6
cd CHEM501_Workshop6
pip install -r requirements.txt
```


