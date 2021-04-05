## referece: https://onehwengineer.com/klayout-python-tutorial-part1/ ##


import pya
import cell_math as cm


class euler_bend_90:

	def __init__(self, num_pts  = 64, layer = 1, cell_name = 'euler_90', radius = 10, angle = 90, p_parameter = 1, use_eff = False, width = 1):
		self.cell_name = cell_name
		self.layer = layer
		self.radius = radius
		self.angle = angle
		self.p_parameter = p_parameter
		self.use_eff = use_eff
		self.num_pts = num_pts
		self.width = width
		self.cell_name = cell_name + '_radius_'+ str(self.radius)

		self.dbu = 1000
		self.layout = pya.Layout()

	def set_p_parameter(self, p):
		self.p_parameter = p

	def set_radius(self, radius):
		if radius < 10:
			raise TypeError('minimum radius is 10um')
		else:
			self.radius = radius

	def set_layer(self, layer):
		if layer < 0 and isinstance(layer, int) == False:
			raise TypeError('layer is an integer larger than 0')
		else:
			self.layer = layer

	def set_num_pts(self, num_pts):
		if num_pts < 4 and isinstance(layer, int) == False:
			raise TypeError('number of points needs to be an integer larger than 4')
		else:
			self.layer = layer

	def set_cell_name(self, cell_name):
			self.cell_name = str(cell_name)

	def generate_shape(self):

		self.UNIT = self.layout.create_cell(str(self.cell_name))

		self.layer = self.layout.layer(self.layer, 0) # Metal


		euler_pts = cm.euler(self.radius, self.angle, self.p_parameter, self.use_eff, self.num_pts)
		pts = []
		for i in range(len(euler_pts)):

			pts.append(pya.Point.from_dpoint(pya.DPoint(euler_pts[i][0]* self.dbu ,euler_pts[i][1]* self.dbu)))

		self.UNIT.shapes(self.layer).insert(pya.Path(pts,self.width*self.dbu))

	def generate_gds(self, name):

		self.UNIT.write(str(name)+ '.gds')



euler = euler_bend_90()
euler.set_radius(20)
euler.set_layer(2)
euler.generate_shape()
euler.generate_gds('euler10')


euler = euler_bend_90()
euler.set_radius(20)
euler.generate_shape()
euler.generate_gds('euler20')


# layout = pya.Layout()


# # Create Cell obj
# UNIT = layout.create_cell("90_euler_bend")


# # Create layer #'s
# l_1x1_outline = layout.layer(1, 0) # 1x1 Outline
# l_metal = layout.layer(11, 0) # Metal


# # Metal dimensions
# line_width = 5*1000 # 10 um
# pitch = 100*1000 # 100 um


# # Draw outline
# # outline = UNIT.shapes(l_1x1_outline).insert( pya.Box(0, 0, pitch, pitch) ) 


# # # Draw metal legs
# # leg1 = UNIT.shapes(l_metal).insert( pya.Box(0, 0, line_width, pitch) ) 
# # leg2 = UNIT.shapes(l_metal).insert( pya.Box(0, pitch-line_width, pitch, pitch) ) 


#     # create the shape

# euler_bend = euler.euler()
# width = 2

# pts = []
# for i in range(len(euler_bend)):

#     pts.append(pya.Point.from_dpoint(pya.DPoint(euler_bend[i][0]* 1000.,euler_bend[i][1]* 1000.)))
# #        pts.append(pya.Point.from_dpoint(pya.DPoint(euler_points[i][0]/self.layout.dbu , euler_points[i][1]/self.layout.dbu )))

# # pt1 = pya.Point.from_dpoint(pya.DPoint(0*1000,0*1000))
# # pt2 = pya.Point.from_dpoint(pya.DPoint(1*1000,1*1000))
# # pt3 = pya.Point.from_dpoint(pya.DPoint(2*1000,2*1000))
# # pt4 = pya.Point.from_dpoint(pya.DPoint(3*1000,3*1000))
# # pts = [pt1,pt2,pt3,pt4]
# leg3 = UNIT.shapes(l_metal).insert(pya.Path(pts,width*1000))




# # Export GDS
# layout.write("0_unit_1x1.gds")
