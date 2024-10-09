#Check Element mesh quality prior to anything else. Then either manually create a set called Faces2 with the elements to examine
# or setup a pressure load and pull from that list and set below. 

#Create set
*createentity sets cardimage=SET_ELEM includeid=0 name="Faces"
*setvalue sets id=3 ids={elems }


#########################

# Code snippet that sets elem_list with all elements in Faces2 set, need to add all elements desired into that 
# group but otherwise is plug and play, can get from pressure panel on whole top of gel, just set face angle = 0 

hm_createmark elems 1 "by set name" "Faces"
set elem_list [hm_getmark elems 1]


#######################

# Section that will be used to turn off the GUI elements so the loop runs much faster. This will obviously break things so will need
# turn it back off when done running a section of code. Also need to call this process. 
namespace eval Panel {


}

proc ::Panel::Performance {flag} {
   switch $flag {
      "on" {
         *retainmarkselections 0
         *entityhighlighting 0
         hm_blockredraw 1
         hm_blockmessages 1
         hm_commandfilestate 0
hwbrowsermanager view flush false
hmbr_signals buffer start
hm_setmouse 0
      }
      "off" -
      default {
         *retainmarkselections 1
         *entityhighlighting 1
         hm_blockredraw 0
         hm_blockmessages 0
         hm_commandfilestate 1
hwbrowsermanager view flush true
hmbr_signals buffer end
hm_setmouse 1
      }
   }
}




#######################

# Simple version if elem_list is already defined THIS WORKS Trick is need to figure out what set of nodes defines 
# the surface and only the surface of the elements in question. Can get from a pressure load on all surface elements. 
# If we use 2D elements (probably the way to go) then we don't need to define nodes just set a mark. That should work. 

::Panel::Performance "on"

set i 1
foreach j $elem_list {

	#puts "$j"
	
	*createmark elems 1 "$j"
	*createmark nodes 1
	*createentity loadcols includeid=0 name="LoadX${j}"
	*pressuresonentity_curve elements 1 1 1 0 0 1 0 1 0 0 0 0 0
	
	*createmark elems 1 "$j"
	*createmark nodes 1
	*createentity loadcols includeid=0 name="LoadY${j}"
	*pressuresonentity_curve elements 1 1 0 1 0 1 0 1 0 0 0 0 0
	
	*createmark elems 1 "$j"
	*createmark nodes 1
	*createentity loadcols includeid=0 name="LoadZ${j}"
	*pressuresonentity_curve elements 1 1 0 0 1 1 0 1 0 0 0 0 0
	
	incr i

}


::Panel::Performance "off"

################# 

#Need to manually set SPCs for the bottom of the gel and change SPCID to that load collector, it will be the last one. 
#Will also need to update loop final number. Then you should be able to run optistruct. 

# We could automate this by creating a group of those nodes and setting them here in a new collector. 

set tempid 1

::Panel::Performance "on"

for { set i 1 } { $i < 8935 } { incr i } {
	


*createentity loadsteps includeid=0 name="Load${i}"

#puts "$tempid"
#puts "$tempid"
#puts "$i"

set spcid 8935
set loadid [expr {$i} ]


################
# Changes the Analysis Type to Linear Static
################
*setvalue loadsteps id=$tempid STATUS=2 OS_TYPE=1
*setvalue loadsteps id=$tempid STATUS=1 4709=1
*setvalue loadsteps id=$tempid STATUS=2 4059=1
*setvalue loadsteps id=$tempid STATUS=2 4060=STATICS

*setvalue loadsteps id=$tempid STATUS=2 3451=0
*setvalue loadsteps id=$tempid STATUS=2 4152=0

################
# Sets the SPCs
################
*setvalue loadsteps id=$tempid STATUS=2 OS_SPCID={loadcols $spcid}
*setvalue loadsteps id=$tempid STATUS=2 4143=1
*setvalue loadsteps id=$tempid STATUS=1 4144=1
*setvalue loadsteps id=$tempid STATUS=1 4145={Loadcols $spcid}
*setvalue loadsteps id=$tempid STATUS=2 ids={$spcid}
*setvalue loadsteps id=1 STATUS=2 ids={1 43 45-82}

################
# Sets the Load
################
*setvalue loadsteps id=$tempid STATUS=2 OS_LOADID={loadcols $loadid}
*setvalue loadsteps id=$tempid STATUS=2 4143=1
*setvalue loadsteps id=$tempid STATUS=1 4146=1
*setvalue loadsteps id=$tempid STATUS=1 4147={Loadcols $loadid}
*setvalue loadsteps id=$tempid STATUS=2 ids={$loadid $spcid}
*setvalue loadsteps id=$tempid STATUS=0 7763=0
*setvalue loadsteps id=$tempid STATUS=1 7740={Loadcols 0}
*setvalue loadsteps id=$tempid STATUS=2 ids={$loadid $spcid}


incr tempid


}

::Panel::Performance "off"