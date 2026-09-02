#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#   Copyright (c) 2002-2022 The ROMS/TOMS Group            Chuning Wang :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::: Benjamin K. Galton-Fenzi :::

local_sub  := ROMS/IceShelf

local_src  := $(wildcard $(local_sub)/*.F)

sources    += $(local_src)

$(eval $(compile-rules))
