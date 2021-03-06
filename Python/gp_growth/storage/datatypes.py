from tables import *

METADATA_TYPES = ["str","int","float"]
PLATE_TYPES = ["Bioscreen","Plate96Well","Plate384Well"]

class Bioscreen(IsDescription):
	time = Float64Col()

	# define a column for each well
	well0 = Float64Col()
	well1 = Float64Col()
	well2 = Float64Col()
	well3 = Float64Col()
	well4 = Float64Col()
	well5 = Float64Col()
	well6 = Float64Col()
	well7 = Float64Col()
	well8 = Float64Col()
	well9 = Float64Col()
	well10 = Float64Col()
	well11 = Float64Col()
	well12 = Float64Col()
	well13 = Float64Col()
	well14 = Float64Col()
	well15 = Float64Col()
	well16 = Float64Col()
	well17 = Float64Col()
	well18 = Float64Col()
	well19 = Float64Col()
	well20 = Float64Col()
	well21 = Float64Col()
	well22 = Float64Col()
	well23 = Float64Col()
	well24 = Float64Col()
	well25 = Float64Col()
	well26 = Float64Col()
	well27 = Float64Col()
	well28 = Float64Col()
	well29 = Float64Col()
	well30 = Float64Col()
	well31 = Float64Col()
	well32 = Float64Col()
	well33 = Float64Col()
	well34 = Float64Col()
	well35 = Float64Col()
	well36 = Float64Col()
	well37 = Float64Col()
	well38 = Float64Col()
	well39 = Float64Col()
	well40 = Float64Col()
	well41 = Float64Col()
	well42 = Float64Col()
	well43 = Float64Col()
	well44 = Float64Col()
	well45 = Float64Col()
	well46 = Float64Col()
	well47 = Float64Col()
	well48 = Float64Col()
	well49 = Float64Col()
	well50 = Float64Col()
	well51 = Float64Col()
	well52 = Float64Col()
	well53 = Float64Col()
	well54 = Float64Col()
	well55 = Float64Col()
	well56 = Float64Col()
	well57 = Float64Col()
	well58 = Float64Col()
	well59 = Float64Col()
	well60 = Float64Col()
	well61 = Float64Col()
	well62 = Float64Col()
	well63 = Float64Col()
	well64 = Float64Col()
	well65 = Float64Col()
	well66 = Float64Col()
	well67 = Float64Col()
	well68 = Float64Col()
	well69 = Float64Col()
	well70 = Float64Col()
	well71 = Float64Col()
	well72 = Float64Col()
	well73 = Float64Col()
	well74 = Float64Col()
	well75 = Float64Col()
	well76 = Float64Col()
	well77 = Float64Col()
	well78 = Float64Col()
	well79 = Float64Col()
	well80 = Float64Col()
	well81 = Float64Col()
	well82 = Float64Col()
	well83 = Float64Col()
	well84 = Float64Col()
	well85 = Float64Col()
	well86 = Float64Col()
	well87 = Float64Col()
	well88 = Float64Col()
	well89 = Float64Col()
	well90 = Float64Col()
	well91 = Float64Col()
	well92 = Float64Col()
	well93 = Float64Col()
	well94 = Float64Col()
	well95 = Float64Col()
	well96 = Float64Col()
	well97 = Float64Col()
	well98 = Float64Col()
	well99 = Float64Col()
	well100 = Float64Col()
	well101 = Float64Col()
	well102 = Float64Col()
	well103 = Float64Col()
	well104 = Float64Col()
	well105 = Float64Col()
	well106 = Float64Col()
	well107 = Float64Col()
	well108 = Float64Col()
	well109 = Float64Col()
	well110 = Float64Col()
	well111 = Float64Col()
	well112 = Float64Col()
	well113 = Float64Col()
	well114 = Float64Col()
	well115 = Float64Col()
	well116 = Float64Col()
	well117 = Float64Col()
	well118 = Float64Col()
	well119 = Float64Col()
	well120 = Float64Col()
	well121 = Float64Col()
	well122 = Float64Col()
	well123 = Float64Col()
	well124 = Float64Col()
	well125 = Float64Col()
	well126 = Float64Col()
	well127 = Float64Col()
	well128 = Float64Col()
	well129 = Float64Col()
	well130 = Float64Col()
	well131 = Float64Col()
	well132 = Float64Col()
	well133 = Float64Col()
	well134 = Float64Col()
	well135 = Float64Col()
	well136 = Float64Col()
	well137 = Float64Col()
	well138 = Float64Col()
	well139 = Float64Col()
	well140 = Float64Col()
	well141 = Float64Col()
	well142 = Float64Col()
	well143 = Float64Col()
	well144 = Float64Col()
	well145 = Float64Col()
	well146 = Float64Col()
	well147 = Float64Col()
	well148 = Float64Col()
	well149 = Float64Col()
	well150 = Float64Col()
	well151 = Float64Col()
	well152 = Float64Col()
	well153 = Float64Col()
	well154 = Float64Col()
	well155 = Float64Col()
	well156 = Float64Col()
	well157 = Float64Col()
	well158 = Float64Col()
	well159 = Float64Col()
	well160 = Float64Col()
	well161 = Float64Col()
	well162 = Float64Col()
	well163 = Float64Col()
	well164 = Float64Col()
	well165 = Float64Col()
	well166 = Float64Col()
	well167 = Float64Col()
	well168 = Float64Col()
	well169 = Float64Col()
	well170 = Float64Col()
	well171 = Float64Col()
	well172 = Float64Col()
	well173 = Float64Col()
	well174 = Float64Col()
	well175 = Float64Col()
	well176 = Float64Col()
	well177 = Float64Col()
	well178 = Float64Col()
	well179 = Float64Col()
	well180 = Float64Col()
	well181 = Float64Col()
	well182 = Float64Col()
	well183 = Float64Col()
	well184 = Float64Col()
	well185 = Float64Col()
	well186 = Float64Col()
	well187 = Float64Col()
	well188 = Float64Col()
	well189 = Float64Col()
	well190 = Float64Col()
	well191 = Float64Col()
	well192 = Float64Col()
	well193 = Float64Col()
	well194 = Float64Col()
	well195 = Float64Col()
	well196 = Float64Col()
	well197 = Float64Col()
	well198 = Float64Col()
	well199 = Float64Col()

class Plate96Well(IsDescription):
	time = Float64Col()

	# define a column for each well
	well0 = Float64Col()
	well1 = Float64Col()
	well2 = Float64Col()
	well3 = Float64Col()
	well4 = Float64Col()
	well5 = Float64Col()
	well6 = Float64Col()
	well7 = Float64Col()
	well8 = Float64Col()
	well9 = Float64Col()
	well10 = Float64Col()
	well11 = Float64Col()
	well12 = Float64Col()
	well13 = Float64Col()
	well14 = Float64Col()
	well15 = Float64Col()
	well16 = Float64Col()
	well17 = Float64Col()
	well18 = Float64Col()
	well19 = Float64Col()
	well20 = Float64Col()
	well21 = Float64Col()
	well22 = Float64Col()
	well23 = Float64Col()
	well24 = Float64Col()
	well25 = Float64Col()
	well26 = Float64Col()
	well27 = Float64Col()
	well28 = Float64Col()
	well29 = Float64Col()
	well30 = Float64Col()
	well31 = Float64Col()
	well32 = Float64Col()
	well33 = Float64Col()
	well34 = Float64Col()
	well35 = Float64Col()
	well36 = Float64Col()
	well37 = Float64Col()
	well38 = Float64Col()
	well39 = Float64Col()
	well40 = Float64Col()
	well41 = Float64Col()
	well42 = Float64Col()
	well43 = Float64Col()
	well44 = Float64Col()
	well45 = Float64Col()
	well46 = Float64Col()
	well47 = Float64Col()
	well48 = Float64Col()
	well49 = Float64Col()
	well50 = Float64Col()
	well51 = Float64Col()
	well52 = Float64Col()
	well53 = Float64Col()
	well54 = Float64Col()
	well55 = Float64Col()
	well56 = Float64Col()
	well57 = Float64Col()
	well58 = Float64Col()
	well59 = Float64Col()
	well60 = Float64Col()
	well61 = Float64Col()
	well62 = Float64Col()
	well63 = Float64Col()
	well64 = Float64Col()
	well65 = Float64Col()
	well66 = Float64Col()
	well67 = Float64Col()
	well68 = Float64Col()
	well69 = Float64Col()
	well70 = Float64Col()
	well71 = Float64Col()
	well72 = Float64Col()
	well73 = Float64Col()
	well74 = Float64Col()
	well75 = Float64Col()
	well76 = Float64Col()
	well77 = Float64Col()
	well78 = Float64Col()
	well79 = Float64Col()
	well80 = Float64Col()
	well81 = Float64Col()
	well82 = Float64Col()
	well83 = Float64Col()
	well84 = Float64Col()
	well85 = Float64Col()
	well86 = Float64Col()
	well87 = Float64Col()
	well88 = Float64Col()
	well89 = Float64Col()
	well90 = Float64Col()
	well91 = Float64Col()
	well92 = Float64Col()
	well93 = Float64Col()
	well94 = Float64Col()
	well95 = Float64Col()

class Plate384Well(IsDescription):
	time = Float64Col()

	# define a column for each well
	well0 = Float64Col()
	well1 = Float64Col()
	well2 = Float64Col()
	well3 = Float64Col()
	well4 = Float64Col()
	well5 = Float64Col()
	well6 = Float64Col()
	well7 = Float64Col()
	well8 = Float64Col()
	well9 = Float64Col()
	well10 = Float64Col()
	well11 = Float64Col()
	well12 = Float64Col()
	well13 = Float64Col()
	well14 = Float64Col()
	well15 = Float64Col()
	well16 = Float64Col()
	well17 = Float64Col()
	well18 = Float64Col()
	well19 = Float64Col()
	well20 = Float64Col()
	well21 = Float64Col()
	well22 = Float64Col()
	well23 = Float64Col()
	well24 = Float64Col()
	well25 = Float64Col()
	well26 = Float64Col()
	well27 = Float64Col()
	well28 = Float64Col()
	well29 = Float64Col()
	well30 = Float64Col()
	well31 = Float64Col()
	well32 = Float64Col()
	well33 = Float64Col()
	well34 = Float64Col()
	well35 = Float64Col()
	well36 = Float64Col()
	well37 = Float64Col()
	well38 = Float64Col()
	well39 = Float64Col()
	well40 = Float64Col()
	well41 = Float64Col()
	well42 = Float64Col()
	well43 = Float64Col()
	well44 = Float64Col()
	well45 = Float64Col()
	well46 = Float64Col()
	well47 = Float64Col()
	well48 = Float64Col()
	well49 = Float64Col()
	well50 = Float64Col()
	well51 = Float64Col()
	well52 = Float64Col()
	well53 = Float64Col()
	well54 = Float64Col()
	well55 = Float64Col()
	well56 = Float64Col()
	well57 = Float64Col()
	well58 = Float64Col()
	well59 = Float64Col()
	well60 = Float64Col()
	well61 = Float64Col()
	well62 = Float64Col()
	well63 = Float64Col()
	well64 = Float64Col()
	well65 = Float64Col()
	well66 = Float64Col()
	well67 = Float64Col()
	well68 = Float64Col()
	well69 = Float64Col()
	well70 = Float64Col()
	well71 = Float64Col()
	well72 = Float64Col()
	well73 = Float64Col()
	well74 = Float64Col()
	well75 = Float64Col()
	well76 = Float64Col()
	well77 = Float64Col()
	well78 = Float64Col()
	well79 = Float64Col()
	well80 = Float64Col()
	well81 = Float64Col()
	well82 = Float64Col()
	well83 = Float64Col()
	well84 = Float64Col()
	well85 = Float64Col()
	well86 = Float64Col()
	well87 = Float64Col()
	well88 = Float64Col()
	well89 = Float64Col()
	well90 = Float64Col()
	well91 = Float64Col()
	well92 = Float64Col()
	well93 = Float64Col()
	well94 = Float64Col()
	well95 = Float64Col()
	well96 = Float64Col()
	well97 = Float64Col()
	well98 = Float64Col()
	well99 = Float64Col()
	well100 = Float64Col()
	well101 = Float64Col()
	well102 = Float64Col()
	well103 = Float64Col()
	well104 = Float64Col()
	well105 = Float64Col()
	well106 = Float64Col()
	well107 = Float64Col()
	well108 = Float64Col()
	well109 = Float64Col()
	well110 = Float64Col()
	well111 = Float64Col()
	well112 = Float64Col()
	well113 = Float64Col()
	well114 = Float64Col()
	well115 = Float64Col()
	well116 = Float64Col()
	well117 = Float64Col()
	well118 = Float64Col()
	well119 = Float64Col()
	well120 = Float64Col()
	well121 = Float64Col()
	well122 = Float64Col()
	well123 = Float64Col()
	well124 = Float64Col()
	well125 = Float64Col()
	well126 = Float64Col()
	well127 = Float64Col()
	well128 = Float64Col()
	well129 = Float64Col()
	well130 = Float64Col()
	well131 = Float64Col()
	well132 = Float64Col()
	well133 = Float64Col()
	well134 = Float64Col()
	well135 = Float64Col()
	well136 = Float64Col()
	well137 = Float64Col()
	well138 = Float64Col()
	well139 = Float64Col()
	well140 = Float64Col()
	well141 = Float64Col()
	well142 = Float64Col()
	well143 = Float64Col()
	well144 = Float64Col()
	well145 = Float64Col()
	well146 = Float64Col()
	well147 = Float64Col()
	well148 = Float64Col()
	well149 = Float64Col()
	well150 = Float64Col()
	well151 = Float64Col()
	well152 = Float64Col()
	well153 = Float64Col()
	well154 = Float64Col()
	well155 = Float64Col()
	well156 = Float64Col()
	well157 = Float64Col()
	well158 = Float64Col()
	well159 = Float64Col()
	well160 = Float64Col()
	well161 = Float64Col()
	well162 = Float64Col()
	well163 = Float64Col()
	well164 = Float64Col()
	well165 = Float64Col()
	well166 = Float64Col()
	well167 = Float64Col()
	well168 = Float64Col()
	well169 = Float64Col()
	well170 = Float64Col()
	well171 = Float64Col()
	well172 = Float64Col()
	well173 = Float64Col()
	well174 = Float64Col()
	well175 = Float64Col()
	well176 = Float64Col()
	well177 = Float64Col()
	well178 = Float64Col()
	well179 = Float64Col()
	well180 = Float64Col()
	well181 = Float64Col()
	well182 = Float64Col()
	well183 = Float64Col()
	well184 = Float64Col()
	well185 = Float64Col()
	well186 = Float64Col()
	well187 = Float64Col()
	well188 = Float64Col()
	well189 = Float64Col()
	well190 = Float64Col()
	well191 = Float64Col()
	well192 = Float64Col()
	well193 = Float64Col()
	well194 = Float64Col()
	well195 = Float64Col()
	well196 = Float64Col()
	well197 = Float64Col()
	well198 = Float64Col()
	well199 = Float64Col()
	well200 = Float64Col()
	well201 = Float64Col()
	well202 = Float64Col()
	well203 = Float64Col()
	well204 = Float64Col()
	well205 = Float64Col()
	well206 = Float64Col()
	well207 = Float64Col()
	well208 = Float64Col()
	well209 = Float64Col()
	well210 = Float64Col()
	well211 = Float64Col()
	well212 = Float64Col()
	well213 = Float64Col()
	well214 = Float64Col()
	well215 = Float64Col()
	well216 = Float64Col()
	well217 = Float64Col()
	well218 = Float64Col()
	well219 = Float64Col()
	well220 = Float64Col()
	well221 = Float64Col()
	well222 = Float64Col()
	well223 = Float64Col()
	well224 = Float64Col()
	well225 = Float64Col()
	well226 = Float64Col()
	well227 = Float64Col()
	well228 = Float64Col()
	well229 = Float64Col()
	well230 = Float64Col()
	well231 = Float64Col()
	well232 = Float64Col()
	well233 = Float64Col()
	well234 = Float64Col()
	well235 = Float64Col()
	well236 = Float64Col()
	well237 = Float64Col()
	well238 = Float64Col()
	well239 = Float64Col()
	well240 = Float64Col()
	well241 = Float64Col()
	well242 = Float64Col()
	well243 = Float64Col()
	well244 = Float64Col()
	well245 = Float64Col()
	well246 = Float64Col()
	well247 = Float64Col()
	well248 = Float64Col()
	well249 = Float64Col()
	well250 = Float64Col()
	well251 = Float64Col()
	well252 = Float64Col()
	well253 = Float64Col()
	well254 = Float64Col()
	well255 = Float64Col()
	well256 = Float64Col()
	well257 = Float64Col()
	well258 = Float64Col()
	well259 = Float64Col()
	well260 = Float64Col()
	well261 = Float64Col()
	well262 = Float64Col()
	well263 = Float64Col()
	well264 = Float64Col()
	well265 = Float64Col()
	well266 = Float64Col()
	well267 = Float64Col()
	well268 = Float64Col()
	well269 = Float64Col()
	well270 = Float64Col()
	well271 = Float64Col()
	well272 = Float64Col()
	well273 = Float64Col()
	well274 = Float64Col()
	well275 = Float64Col()
	well276 = Float64Col()
	well277 = Float64Col()
	well278 = Float64Col()
	well279 = Float64Col()
	well280 = Float64Col()
	well281 = Float64Col()
	well282 = Float64Col()
	well283 = Float64Col()
	well284 = Float64Col()
	well285 = Float64Col()
	well286 = Float64Col()
	well287 = Float64Col()
	well288 = Float64Col()
	well289 = Float64Col()
	well290 = Float64Col()
	well291 = Float64Col()
	well292 = Float64Col()
	well293 = Float64Col()
	well294 = Float64Col()
	well295 = Float64Col()
	well296 = Float64Col()
	well297 = Float64Col()
	well298 = Float64Col()
	well299 = Float64Col()
	well300 = Float64Col()
	well301 = Float64Col()
	well302 = Float64Col()
	well303 = Float64Col()
	well304 = Float64Col()
	well305 = Float64Col()
	well306 = Float64Col()
	well307 = Float64Col()
	well308 = Float64Col()
	well309 = Float64Col()
	well310 = Float64Col()
	well311 = Float64Col()
	well312 = Float64Col()
	well313 = Float64Col()
	well314 = Float64Col()
	well315 = Float64Col()
	well316 = Float64Col()
	well317 = Float64Col()
	well318 = Float64Col()
	well319 = Float64Col()
	well320 = Float64Col()
	well321 = Float64Col()
	well322 = Float64Col()
	well323 = Float64Col()
	well324 = Float64Col()
	well325 = Float64Col()
	well326 = Float64Col()
	well327 = Float64Col()
	well328 = Float64Col()
	well329 = Float64Col()
	well330 = Float64Col()
	well331 = Float64Col()
	well332 = Float64Col()
	well333 = Float64Col()
	well334 = Float64Col()
	well335 = Float64Col()
	well336 = Float64Col()
	well337 = Float64Col()
	well338 = Float64Col()
	well339 = Float64Col()
	well340 = Float64Col()
	well341 = Float64Col()
	well342 = Float64Col()
	well343 = Float64Col()
	well344 = Float64Col()
	well345 = Float64Col()
	well346 = Float64Col()
	well347 = Float64Col()
	well348 = Float64Col()
	well349 = Float64Col()
	well350 = Float64Col()
	well351 = Float64Col()
	well352 = Float64Col()
	well353 = Float64Col()
	well354 = Float64Col()
	well355 = Float64Col()
	well356 = Float64Col()
	well357 = Float64Col()
	well358 = Float64Col()
	well359 = Float64Col()
	well360 = Float64Col()
	well361 = Float64Col()
	well362 = Float64Col()
	well363 = Float64Col()
	well364 = Float64Col()
	well365 = Float64Col()
	well366 = Float64Col()
	well367 = Float64Col()
	well368 = Float64Col()
	well369 = Float64Col()
	well370 = Float64Col()
	well371 = Float64Col()
	well372 = Float64Col()
	well373 = Float64Col()
	well374 = Float64Col()
	well375 = Float64Col()
	well376 = Float64Col()
	well377 = Float64Col()
	well378 = Float64Col()
	well379 = Float64Col()
	well380 = Float64Col()
	well381 = Float64Col()
	well382 = Float64Col()
	well383 = Float64Col()

class Metadata(IsDescription):
	name = StringCol(16)
	val = StringCol(16)
	type = EnumCol(METADATA_TYPES, 'str', base='uint8')
	well_array = StringCol(16)