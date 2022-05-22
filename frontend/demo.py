#!/usr/bin/env python3
import tkinter as tk
from LMChips import LMChips

def main():
	width = 800
	height = 800
	#instance of tkinter
	root = tk.Tk()
	root.geometry(f"{width}x{height+50}")
	root.title("Chips Demo")

	#main frame and canvas in frame
	frame = tk.Frame(root)
	canvas = tk.Canvas(frame, width = width, height = height)
	canvas.pack()

	LM = LMChips(canvas, width, height)
	LM.update_triggers_refresh = False
	#defaut shape
	LM.lines = [[
		(.45, .7), (-.99, .0),
		(.15, .9), ( .5,  .5),
		(.45, .1), ( .99, .0),
		(.75, .9), (-.99, .0),
		(.75, .1), ( .99, .0),
	]]
	LM.refresh()

	buttons = {
		"new_line":	tk.Button(frame, text = "START NEW LINE", command = LM.start_new_line),
		"refresh":	tk.Button(frame, text = "REFRESH", command = LM.refresh),
		"undo":		tk.Button(frame, text = "UNDO", command = LM.undo),
		"redo":		tk.Button(frame, text = "REDO", command = LM.redo),
	}
	for button in buttons.values():
		button.pack(side = tk.LEFT)

	# # scrol_y = tk.Scrollbar(frame, orient = tk.VERTICAL)
	# # Adding text widget for inserting images
	# pdf = tk.Text(frame)
	# image_item = pdf.image_create(tk.END, image = image)
	# pdf.pack(fill = tk.BOTH, expand = 1)
	frame.pack()

	root.mainloop()

if __name__ == '__main__':
	main()
