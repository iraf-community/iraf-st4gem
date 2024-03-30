      real function unitrp (nwave, wave, band, units)
      integer nwave
      real wave(*)
      real band(*)
      integer*2 units(*)
      integer utype
      logical Memb(1)
      integer*2 Memc(1)
      integer*2 Mems(1)
      integer Memi(1)
      integer*4 Meml(1)
      real Memr(1)
      double precision Memd(1)
      complex Memx(1)
      equivalence (Memb, Memc, Mems, Memi, Meml, Memr, Memd, Memx)
      common /Mem/ Memd
      integer sp
      integer form
      integer vflux
      real sum
      real resp
      integer wordmh
      real sumfit
      real asumr
      real cntrae
      integer sw0001
      integer*2 baduns(19)
      integer*2 formlt(64)
      logical xerflg
      common /xercom/ xerflg
      integer*2 st0001(8)
      save
      integer iyy
      data (baduns(iyy),iyy= 1, 8) / 85,110,107,110,111,119,110, 32/
      data (baduns(iyy),iyy= 9,16) /102,108,117,120, 32,117,110,105/
      data (baduns(iyy),iyy=17,19) /116,115, 0/
      data (formlt(iyy),iyy= 1, 8) /112,104,111,116,108, 97,109, 44/
      data (formlt(iyy),iyy= 9,16) / 99,111,117,110,116,115, 44,102/
      data (formlt(iyy),iyy=17,24) /108, 97,109, 44,102,110,117, 44/
      data (formlt(iyy),iyy=25,32) /112,104,111,116,110,117, 44,106/
      data (formlt(iyy),iyy=33,40) /121, 44,109,106,121, 44, 97, 98/
      data (formlt(iyy),iyy=41,48) /109, 97,103, 44,115,116,109, 97/
      data (formlt(iyy),iyy=49,56) /103, 44,118,101,103, 97,109, 97/
      data (formlt(iyy),iyy=57,64) /103, 44,111, 98,109, 97,103, 0/
      data st0001 /112,104,111,116,108, 97,109, 0/
      unitrp = 0
         call smark (sp)
         call salloc (form, 63 , 2)
         resp = 0.0
         call xstrcy(units, memc(form), 63 )
         call strfix (memc(form))
         utype = wordmh (memc(form), formlt)
         sw0001=(utype)
         goto 110
120      continue
            sum = sumfit (nwave, wave, 0, band)
            if (.not.(sum .gt. 0.0)) goto 130
               resp = 1.0 / sum
130         continue
         goto 111
140      continue
            sum = asumr (band, nwave)
            if (.not.(sum .gt. 0.0)) goto 150
               resp = 1.0 / sum
150         continue
         goto 111
160      continue
            sum = sumfit (nwave, wave, 1, band)
            sum = sum / (6.62620e-27 * 2.997925e18 )
            if (.not.(sum .gt. 0.0)) goto 170
               resp = 1.0 / sum
170         continue
         goto 111
180      continue
            sum = sumfit (nwave, wave, -1, band)
            sum = sum / 6.62620e-27
            if (.not.(sum .gt. 0.0)) goto 190
               resp = 1.0 / sum
190         continue
         goto 111
200      continue
            sum = sumfit (nwave, wave, -2, band)
            sum = sum * 2.997925e18
            if (.not.(sum .gt. 0.0)) goto 210
               resp = 1.0 / sum
210         continue
         goto 111
220      continue
            sum = sumfit (nwave, wave, -1, band)
            sum = sum * 1.0e-23 / 6.62620e-27
            if (.not.(sum .gt. 0.0)) goto 230
               resp = 1.0 / sum
230         continue
         goto 111
240      continue
            sum = sumfit (nwave, wave, -1, band)
            sum = sum * 1.0e-26 / 6.62620e-27
            if (.not.(sum .gt. 0.0)) goto 250
               resp = 1.0 / sum
250         continue
         goto 111
260      continue
            sum = sumfit (nwave, wave, -1, band)
            sum = sum / 6.62620e-27
            if (.not.(sum .gt. 0.0)) goto 270
               resp = 2.5 * alog10 (sum) + (-48.60)
270         continue
         goto 111
280      continue
            sum = sumfit (nwave, wave, 1, band)
            sum = sum / (6.62620e-27 * 2.997925e18 )
            if (.not.(sum .gt. 0.0)) goto 290
               resp = 2.5 * alog10 (sum) + (-21.10)
290         continue
         goto 111
300      continue
            call salloc (vflux, nwave, 6)
            call getvea (nwave, wave, memr(vflux))
            if (xerflg) goto 100
            sum = cntrae (nwave, wave, band, memr(vflux), st0001)
            if (.not.(sum .gt. 0.0)) goto 310
               resp = 2.5 * alog10 (sum)
310         continue
         goto 111
320      continue
            sum = asumr (band, nwave)
            if (.not.(sum .gt. 0.0)) goto 330
               resp = 2.5 * alog10 (sum)
330         continue
         goto 111
340      continue
            call synphr (baduns, units)
            if (xerflg) goto 100
            goto 111
110      continue
            if (sw0001.lt.1.or.sw0001.gt.11) goto 340
            goto (120,140,160,180,200,220,240,260,280,300,320),sw0001
111      continue
         call sfree (sp)
         unitrp = (resp)
         goto 100
100      return
      end
c     synphr  synphoterr
c     getvea  getvega
c     formlt  formlist
c     cntrae  cntrate
c     unitrp  unitresp
c     wordmh  word_match
c     sumfit  sumfilt
c     baduns  badunits
