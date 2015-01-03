from asynch_py.asynch_interface import asynchsolver
import sys

#Parse command line arguments
numargs = len(sys.argv)
if numargs != 2:
	print 'Need an input .gbl file'
	sys.exit(1)

#Prepare system
asynch = asynchsolver()
asynch.Parse_GBL(sys.argv[1])
asynch.Load_System()
N = asynch.Get_Number_Links()
print 'I see',N,'links.'

#Prepare outputs
asynch.Prepare_Temp_Files()
asynch.Write_Current_Step()
asynch.Prepare_Peakflow_Output()
asynch.Prepare_Output()

#Advance solver
asynch.Advance(True)
print 'Calculations done!'

#Take a snapshot
asynch.Take_System_Snapshot(None)

#Create output files
asynch.Create_Output(None)
asynch.Create_Peakflows_Output()

#Cleanup
#asynch.Delete_Temporary_Files()


