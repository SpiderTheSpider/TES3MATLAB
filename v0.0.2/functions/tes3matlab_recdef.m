function rec_types = tes3matlab_recdef
% Defines record types in a TES3 file (besides HEDR and MAST, which are
% special cases)

    rec_types = {'ACTI',...     % ENTITY,   ACTIVATOR
                 'ALCH',...     % ITEM,     POTION / ALCHEMY PRODUCT
                 'APPA',...     % ITEM,     ALCHEMY APPARATUS
                 'ARMO',...     % ITEM,     ARMOR
                 'BODY',...     % (N)PC,    BODY PART (FOR HUMANOIDS)
                 'BOOK',...     % ITEM,     BOOK
                 'BSGN',...     % (N)PC,    BIRTHSIGN
                 'CELL',...     % WORLD,    CELL, INTERIOR OR EXTERIOR
                 'CLAS',...     % (N)PC,    CLASS
                 'CLOT',...     % ITEM,     CLOTHING
                 'CONT',...     % ENTITY,   CONTAINER
                 'CREA',...     % ENTITY,   CREATURE
                 'DIAL',...     % DIALOGUE, TOPIC
                 'DOOR',...     % ENTITY,   DOOR
                 'ENCH',...     % MAGIC,    ENCHANTMENT
                 'FACT',...     % (N)PC,    FACTION
                 'GLOB',...     % VARIABLE, GLOBAL
                 'GMST',...     % VARIABLE, GAME SETTING
                 'INFO',...     % DIALOGUE, RESPONSE
                 'INGR',...     % ITEM,     INGREDIENT
                 'LAND',...     % WORLD,    LANDSCAPE
                 'LEVC',...     % ENTITY,   LEVELED CREATURE
                 'LEVI',...     % ITEM,     LEVELED ITEM
                 'LIGH',...     % ENTITY,   LIGHT
                 'LOCK',...     % ITEM,     LOCKPICKING TOOL
                 'LTEX',...     % WORLD,    LANDSCAPE TEXTURE
                 'MGEF',...     % MAGIC,    MAGIC EFFECT
                 'MISC',...     % ITEM,     MISCELLANEOUS
                 'NPC_',...     % ENTITY,   NPC
                 'PGRD',...     % WORLD,    PATHGRID
                 'PROB',...     % ITEM,     TRAP PROBE
                 'RACE',...     % (N)PC,    RACE
                 'REGN',...     % WORLD,    REGION
                 'REPA',...     % ITEM,     REPAIR TOOL
                 'SCPT',...     % OTHER,    SCRIPT
                 'SKIL',...     % (N)PC,    SKILL
                 'SNDG',...     % ENTITY,   SOUND GENERATOR
                 'SOUN',...     % OTHER,    SOUND
                 'SPEL',...     % MAGIC,    SPELL
                 'SSCR',...     % OTHER,    START SCRIPT
                 'STAT',...     % ENTITY,	STATIC
                 'WEAP',...     % ITEM,     WEAPON
                };

end