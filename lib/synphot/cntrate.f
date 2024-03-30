      real function cntrae (nwave, wave, band, spec, units)
      integer nwave
      real wave(*)
      real band(*)
      real spec(*)
      integer*2 units(*)
      integer utype
      integer done
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
      integer prod
      integer form
      real rate
      integer ismagt
      integer wordmh
      integer anytot
      integer photty
      real sumfit
      real asumr
      integer sw0001
      integer*2 baduns(19)
      integer*2 formlt(64)
      logical xerflg
      common /xercom/ xerflg
      integer*2 st0001(4)
      integer*2 st0002(5)
      integer*2 st0003(5)
      integer*2 st0004(7)
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
      data st0001 /102,110,117, 0/
      data st0002 /102,108, 97,109, 0/
      data st0003 /102,108, 97,109, 0/
      data st0004 / 99,111,117,110,116,115, 0/
      cntrae = 0
         call smark (sp)
         call salloc (prod, nwave, 6)
         call salloc (form, 63 , 2)
         if (.not.(ismagt (units) .eq. 0)) goto 110
            call amulr (band, spec, memr(prod), nwave)
            goto 111
110      continue
            call amovr (spec, memr(prod), nwave)
            done = anytot (units, nwave, wave, memr(prod))
            if (xerflg) goto 100
            call amulr (band, memr(prod), memr(prod), nwave)
111      continue
         call xstrcy(units, memc(form), 63 )
         call strfix (memc(form))
         utype = wordmh (memc(form), formlt)
         sw0001=(utype)
         goto 120
130      continue
            rate = sumfit (nwave, wave, 0, memr(prod))
         goto 121
140      continue
            rate = asumr (memr(prod), nwave)
         goto 121
150      continue
            rate = sumfit (nwave, wave, 1, memr(prod))
            rate = rate / (6.62620e-27 * 2.997925e18 )
         goto 121
160      continue
            rate = sumfit (nwave, wave, -1, memr(prod))
            rate = rate / 6.62620e-27
         goto 121
170      continue
            rate = sumfit (nwave, wave, -2, memr(prod))
            rate = rate * 2.997925e18
         goto 121
180      continue
            rate = sumfit (nwave, wave, -1, memr(prod))
            rate = rate * 1.0e-23 / 6.62620e-27
         goto 121
190      continue
            rate = sumfit (nwave, wave, -1, memr(prod))
            rate = rate * 1.0e-26 / 6.62620e-27
         goto 121
200      continue
            done = photty (st0001, nwave, wave, memr(prod))
            if (xerflg) goto 100
            rate = sumfit (nwave, wave, -1, memr(prod))
            rate = rate / 6.62620e-27
         goto 121
210      continue
            done = photty (st0002, nwave, wave, memr(prod))
            if (xerflg) goto 100
            rate = sumfit (nwave, wave, 1, memr(prod))
            rate = rate / (6.62620e-27 * 2.997925e18 )
         goto 121
220      continue
            done = photty (st0003, nwave, wave, memr(prod))
            if (xerflg) goto 100
            rate = sumfit (nwave, wave, 1, memr(prod))
            rate = rate / (6.62620e-27 * 2.997925e18 )
         goto 121
230      continue
            done = photty (st0004, nwave, wave, memr(prod))
            if (xerflg) goto 100
            rate = asumr (memr(prod), nwave)
         goto 121
240      continue
            call synphr (baduns, units)
            if (xerflg) goto 100
            goto 121
120      continue
            if (sw0001.lt.1.or.sw0001.gt.11) goto 240
            goto (130,140,150,160,170,180,190,200,210,220,230),sw0001
121      continue
         call sfree (sp)
         cntrae = (rate)
         goto 100
100      return
      end
c     synphr  synphoterr
c     anytot  anytophot
c     photty  phottoany
c     formlt  formlist
c     cntrae  cntrate
c     wordmh  word_match
c     ismagt  is_magunit
c     sumfit  sumfilt
c     baduns  badunits
