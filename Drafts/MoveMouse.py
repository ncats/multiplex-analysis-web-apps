import pyautogui
import time

# Set the interval in seconds
# interval = 20

# try:
#     while True:
#         # Move the mouse cursor to a new position
#         pyautogui.moveRel(10, 10, duration=0.5)  # Move 10 pixels right and 10 pixels down
#         time.sleep(interval)
#         pyautogui.moveRel(-10, -10, duration=0.5)  # Move 10 pixels left and 10 pixels up
#         time.sleep(interval)  # Wait for the specified interval
# except KeyboardInterrupt:
#     print("\nProgram terminated.")


# from pynput.mouse import Controller
# import time

# # Set the interval in seconds
# interval = 20

# mouse = Controller()

# try:
#     while True:
#         # Move the mouse cursor to a new position
#         mouse.move(10, 10)
#         time.sleep(interval)  # Wait for the specified interval
# except KeyboardInterrupt:
#     print("\nProgram terminated.")

import pyautogui
import time

while True:
    pyautogui.press('shift')  # Press the shift key
    time.sleep(20)  