def torsion(a, b, c, d):
#    sqrt=math.sqrt, acos=math.acos

    dr=180./pi
    v12x, v12y, v12z = a[0] - b[0], a[1]-b[1], a[2]-b[2]
    v32x, v32y, v32z = c[0] - b[0], c[1]-b[1], c[2]-b[2]
    v43x, v43y, v43z = d[0] - c[0], d[1]-c[1], d[2]-c[2]

    vn13x = v12y*v32z - v12z*v32y
    vn13y = v12z*v32x - v12x*v32z
    vn13z = v12x*v32y - v12y*v32x

    vn24x = v32z*v43y - v32y*v43z
    vn24y = v32x*v43z - v32z*v43x
    vn24z = v32y*v43x - v32x*v43y

    v12 = vn13x*vn24x + vn13y*vn24y + vn13z*vn24z
    v11 = vn13x**2 + vn13y**2 + vn13z**2
    v22 = vn24x**2 + vn24y**2 + vn24z**2
#    print v11*v22,v12,v11,v22
    v1122=v11*v22
    if v1122<1e-16:
        if v12>0:
            return 0.0
        else:
            return -180.0
    else:
        ang = v12/sqrt(v11*v22)
        if ang >= 1.0:
            return 0.0
        elif ang <= -1.0:
            return -180.0
        else:
            ang = acos(ang) * dr

    vtmp = vn13x * (vn24y*v32z - vn24z*v32y) + \
           vn13y * (vn24z*v32x - vn24x*v32z) + \
           vn13z * (vn24x*v32y - vn24y*v32x) < 0.0
    if vtmp:
        return -ang
    else:
        return ang
