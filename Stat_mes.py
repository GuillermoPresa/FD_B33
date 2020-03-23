import cg_pos


BlockFuel = 1000 #[Kg]

class DataBlock:
	
	def __init__(self, PayloadList, PayloadList2 = None):
		self.PayloadList = PayloadList
		self.PayloadList2 = PayloadList2
		self.DataLineList =  []
	

class DataLine:
	
	def __init__(self, hp, IAS, a, de, detr, Fe, FFI, FFr, F_used, TAT, PayloadList):
		self.hp = hp
		self.IAS = IAS
		self.a = a
		self.FFI = FFI
		self.FFr = FFr
		self.FuelUsed = FuelUsed
		self.Fuel = BlockFuel - FuelUsed
		self.TAT = TAT
		self.cg_pos = cg_pos.cg(self.Fuel, PayloadList) 




