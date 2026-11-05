
#include "ArduinoMqttClient.h"
#include "Arduino.h"
#include "Arduino_BHY2Host.h"
#include "WiFiNINA.h"
#include "arduino_secrets.h"

///////please enter your sensitive data in the Secret tab/arduino_secrets.h
char ssid[] = SECRET_SSID;        // your network SSID (name)
char pass[] = SECRET_PASS;    // your network password (use for WPA, or use as key for WEP)

WiFiClient wifiClient;
MqttClient mqttClient(wifiClient);

const char broker[] = SECRET_SERVER;
int        port     = 1883;
const char topic1[] = "joetopic1";
const char topic2[] = "joetopic2";

//set interval for sending messages (milliseconds)
const long interval = 5000;
unsigned long previousMillis = 0;

int count = 0;

void setup() {
  // attempt to connect to Wifi network:
  Serial.print("Attempting to connect to WPA SSID: ");
  Serial.println(ssid);
  while (WiFi.begin(ssid, pass) != WL_CONNECTED) {
    // failed, retry
    delay(5000);
  }

  if (!mqttClient.connect(broker, port)) {
    while (1);
  }

}

void loop() {
  // call poll() regularly to allow the library to send MQTT keep alive which
  // avoids being disconnected by the broker
  mqttClient.poll();

  unsigned long currentMillis = millis();

  if (currentMillis - previousMillis >= interval) {
    // save the last time a message was sent
    previousMillis = currentMillis;

    //record random value from A0, A1 and A2
    int Rvalue1 = analogRead(A0);
    int Rvalue2 = analogRead(A1);

    Serial.print("Sending message to topic1: ");
    Serial.println(topic1);
    Serial.println(Rvalue1);

    Serial.print("Sending message to topic2: ");
    Serial.println(topic2);
    Serial.println(Rvalue2);

    // send message, the Print interface can be used to set the message contents
    mqttClient.beginMessage(topic1);
    mqttClient.print(Rvalue1);
    mqttClient.endMessage();

    mqttClient.beginMessage(topic2);
    mqttClient.print(Rvalue2);
    mqttClient.endMessage();

    Serial.println();
  }
}
