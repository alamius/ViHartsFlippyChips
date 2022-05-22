#!/usr/bin/env python3
import tkinter as tk

class LineManager(object):
	def __init__(self, canvas, width, height, filename = "out"):
		super(LineManager, self).__init__()
		self.canvas = canvas
		self.lines = [[]]
		self.filename = filename
		self.width = width
		self.height = height
		self.clipping_distance = .05
		#this list may contain points (2-tuples) and [] (the latter of which signals to start a new line)
		#it is filled when calling undo and reduced when calling redo, empied when updating with a click even
		self.undo_buffer = []
		self.update_triggers_refresh = True

	def point_near_in_list(self, point, clipping_distance = None):
		"""returns the line index and the point index as `(l, p)`
			if the point `self.lines[l][p]` is within a `clipping_distance` radius of `point`.
		"""

		if clipping_distance is None:
			clipping_distance = self.clipping_distance

		for l in range(len(self.lines)):
			for p in range(len(self.lines[l])):
				other = self.lines[l][p]
				if abs(other[0] - point[0]) > clipping_distance:
					continue

				if abs(other[1] - point[1]) > clipping_distance:
					continue

				if ((other[0] - point[0]) ** 2 + (other[1] - point[1])**2)**.5 < clipping_distance:
					return l, p

	def add_point(self, x, y):
		"""add the point (x, y) to the current line. clipping to a near point is checked."""

		point = (x, y)
		is_near = self.point_near_in_list(point)

		if is_near is None:
			self.lines[-1] += [point]
		else:
			l, p = is_near #unpack the line and point indices that clipping check found
			self.lines[-1] += [self.lines[l][p]]

	def pop_point(self):
		"""pops the last point off the last line.
			if the last line is empty, remove it (doesn't pop off the previous automatically).
		returns the popped value: point aka tuple(number, number), empty list or None
		"""

		if len(self.lines[-1]):
			return self.lines[-1].pop()
		elif len(self.lines):
			return self.lines.pop()
		else:
			return None

	def update(self, event = None):
		"""write the image and bind the click listener.
		if event is given: add its x and y coordinates as a point ot the current list
		"""

		if event is not None:
			self.add_point(event.x / self.width, 1 - event.y / self.height)
			self.undo_buffer = []

		if self.update_triggers_refresh:
			self.refresh()

		return self

	def refresh(self):
		self.write()
		#protect against GC
		self.image = self.get_image()
		self.image_item = self.canvas.create_image((0, 0), anchor = tk.NW, image = self.image)
		self.canvas.tag_bind(self.image_item, '<Button-1>', self.update)

	def start_new_line(self):
		self.lines += [[]]
		self.undo_buffer = []

	def toggle_smoothed(self):
		self.smoothed = not self.smoothed
		self.update()

	def undo(self):
		self.undo_buffer += [self.pop_point()]
		self.update()

	def redo(self):
		if self.undo_buffer == []:
			return

		popped = self.undo_buffer.pop()
		if popped == []:
			self.lines += [[]]
		else:
			if self.lines == []:
				self.lines += [[]]
			self.lines[-1] += [popped]

		self.update()
