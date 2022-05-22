#!/usr/bin/env python3
from PIL import Image, ImageTk
from LineManager import LineManager
from os import system

class LMChips(LineManager):
	def __init__(self, canvas, width, height, filename = "out"):
		super(LMChips, self).__init__(canvas, width, height, filename)

	def add_point(self, x, y):
		"""add the point (x, y) to the current line. clipping to a near point is checked."""

		if len(self.lines[-1]) % 2 == 1:
			X, Y = self.lines[-1][-1]
			point = (3*(x - X), 3*(y - Y))
			self.lines[-1] += [point]
		else:
			super(LMChips, self).add_point(x, y)

	def convert_to_format(self, lines = None):
		"""converts the lists of points that are the lines into pyx paths,
		taking some extra steps to protect against smoothing corruption
		"""

		if lines is None:
			lines = self.lines

		def format(num, lower_limit = 0, lower_limit_numstr = "0"):
			result = ""
			if num < 0:
				result += "-"
			else:
				result += "."

			if num <= lower_limit:
				numstr = lower_limit_numstr
			elif num >= 1:
				numstr = "99"
			else:
				numstr = str(num)
				numstr = numstr[numstr.index(".") + 1 : ]

			result += numstr
			return result

		result = ""
		for points in lines:
			p = 0
			while p < len(points) - 1:
				result += format(points[p][0])
				result += format(points[p][1])
				p += 1
				result += "->"
				result += format(points[p][0], -1, "99")
				result += format(points[p][1], -1, "99")
				p += 1
				result += "\n"

		return result

	def get_image(self):
		system(f"convert {self.filename}.ppm -resize {self.width}x{self.height} {self.filename}.jpg")
		img = Image.open(self.filename + ".jpg")
		# img.resize((self.width, self.height), Image.ANTIALIAS)
		return ImageTk.PhotoImage(img)

	def write(self):
		with open(self.filename + ".chip", "w") as file:
			file.write(self.convert_to_format())

		command = f"../chips {self.filename} --chip-file={self.filename}.chip --color"
		print(command)
		system(command)
