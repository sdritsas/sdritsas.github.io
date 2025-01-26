try:
    ''' Within Grasshopper
    '''
    from Rhino.Geometry import Point3d, BoundingBox, Mesh
    def Show( mesh ):
        pass
except:
    ''' Stand Alone Version
    '''
    import math
    class Point3d:
        def __init__( self, x, y, z ):
            self.X = x
            self.Y = y
            self.Z = z
        def __add__( self, this ):
            return Point3d(
                self.X + this.X,
                self.Y + this.Y,
                self.Z + this.Z )
        def __mul__( self, this ):
            if( isinstance( this, float ) ):
                return Point3d(
                    self.X * this,
                    self.Y * this,
                    self.Z * this )
            return ( self.X * this +
                     self.Y * this +
                     self.Z * this )
        def DistanceTo( self, point ):
            return math.sqrt(
                ( self.X - point.X ) ** 2 +
                ( self.Y - point.Y ) ** 2 +
                ( self.Z - point.Z ) ** 2 )
    class Mesh:
        class ArrayList:
            def __init__( self ):
                self.Array = []
            def Add( self, point ):
                index = len( self.Array )
                self.Array.append( [point.X, point.Y, point.Z] )
                return index
            def AddFace( self, a, b, c, d = None ):
                if( d is None ):
                    self.Array.append( [a, b, c] )
                else:
                    self.Array.append( [a, b, c] )
                    self.Array.append( [c, d, a] )
            @property
            def Count( self ):
                return len( self.Array )
        def __init__( self ):
            self.Vertices = Mesh.ArrayList( )
            self.Faces    = Mesh.ArrayList( )

    Domain, Resolution = 2.0, 64

    def Show( mesh ):
        import trimesh, numpy as np
        mesh = trimesh.Trimesh(
            np.array( mesh.Vertices.Array, dtype=np.float64 ),
            np.array( mesh.Faces.Array, dtype=np.int32 ) )
        mesh.fix_normals( )
        mesh.show( resolution = ( 512, 512 ))

class VoxelGrid:
    def __init__( self, minimum, maximum, resolution ):
        self.Min    = minimum
        self.Max    = maximum
        self.Count  = resolution
        self.Reset( )

    def Reset( self ):
        nx, ny, nz = self.Count
        self.Values = [[[0.0
            for iz in range( nz )]
                for iy in range( ny )]
                    for ix in range( nx )]

    def PointIndices( self ):
        nx, ny, nz = self.Count
        for ix in range( nx ):
            for iy in range( ny ):
                for iz in range( nz ):
                    yield ( ix, iy, iz )

    def ValueAt( self, ix, iy, iz ):
        return self.Values[ix][iy][iz]

    def PointAt( self, ix, iy, iz ):
        nx, ny, nz = self.Count
        tx = ix / ( nx - 1 )
        ty = iy / ( ny - 1 )
        tz = iz / ( nz - 1 )
        return Point3d(
            self.Min.X * ( 1 - tx ) + self.Max.X * tx,
            self.Min.Y * ( 1 - ty ) + self.Max.Y * ty,
            self.Min.Z * ( 1 - tz ) + self.Max.Z * tz )

    def Points( self ):
        for ix, iy, iz in self.PointIndices( ):
            yield self.PointAt( ix, iy, iz )

    def VoxelIndices( self ):
        nx, ny, nz = self.Count
        for ix in range( nx - 1 ):
            for iy in range( ny - 1 ):
                for iz in range( nz - 1 ):
                    yield ( ix, iy, iz )

    def VoxelBoxAt( self, ix, iy, iz ):
        return BoundingBox(
            self.PointAt( ix + 0, iy + 0, iz + 0 ),
            self.PointAt( ix + 1, iy + 1, iz + 1 ) )

    Vertices = [
        [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
        [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]

    def VoxelAt( self, ix, iy, iz ):
        return [self.PointAt( ix + dx, iy + dy, iz + dz )
            for dx, dy, dz in VoxelGrid.Vertices]

    def Voxels( self ):
        for ix, iy, iz in self.VoxelIndices( ):
            yield self.VoxelAt( ix, iy, iz )

    def VoxelBoxAt( self, ix, iy, iz ):
        return BoundingBox(
            self.PointAt( ix + 0, iy + 0, iz + 0 ),
            self.PointAt( ix + 1, iy + 1, iz + 1 ) )

    def VoxelBoxes( self ):
        for ix, iy, iz in self.VoxelIndices( ):
            yield self.VoxelBoxAt( ix, iy, iz )

    def TetrahedronPointsAt( self, ix, iy, iz, index ):
        offsets = self.TetrahedronOffsetsAt( ix, iy, iz, index )
        return [self.PointAt( ix + dx, iy + dy, iz + dz )
                    for dx, dy, dz in offsets]

    def TetrahedronValuesAt( self, ix, iy, iz, index ):
        offsets = self.TetrahedronOffsetsAt( ix, iy, iz, index )
        return [self.Values[ix + dx][iy + dy][iz + dz]
                    for dx, dy, dz in offsets]

    def TetrahedronPointsAndValuesAt( self, ix, iy, iz, index ):
        offsets = self.TetrahedronOffsetsAt( ix, iy, iz, index )
        points, values = [], []
        for dx, dy, dz in offsets:
            points.append( self.PointAt( ix + dx, iy + dy, iz + dz ) )
            values.append( self.Values[ix + dx][iy + dy][iz + dz] )
        return ( points, values )

    Tetras = [
        [[0, 2, 5, 7], [0, 4, 5, 7], [0, 2, 3, 7], [2, 5, 6, 7], [0, 1, 2, 5]],
        [[0, 1, 3, 4], [1, 4, 5, 6], [3, 4, 6, 7], [1, 2, 3, 6], [1, 3, 4, 6]]]

    def TetrahedraAt( self, ix, iy, iz ):
        for index in range( len( VoxelGrid.Tetras[0] ) ):
            yield self.TetrahedronPointsAndValuesAt( ix, iy, iz, index )

    def TetrahedronOffsetsAt( self, ix, iy, iz, index ):
        tetrahedron = self.TetrahedronAt( ix, iy, iz, index )
        return [VoxelGrid.Vertices[vertex]
                    for vertex in tetrahedron]

    def TetrahedronAt( self, ix, iy, iz, index ):
        tetrahedra = self.ReflectionAt( ix, iy, iz )
        return tetrahedra[index]

    def ReflectionAt( self, ix, iy, iz ):
        parity = ( ix ^ iy ^ iz ) & 1
        return VoxelGrid.Tetras[parity]

    TetraFaces = [[0, 1, 2],
                  [1, 2, 3],
                  [2, 0, 3],
                  [0, 1, 3]]

    def TetrahedronMeshAt( self, ix, iy, iz, index ):
        mesh = Mesh( )
        for vertex in  self.TetrahedronPointsAt( ix, iy, iz, index ):
            mesh.Vertices.Add( vertex )
        for a, b, c in VoxelGrid.TetraFaces:
            mesh.Faces.AddFace( a, b, c )
        return mesh

    def VoxelTetrahedraMeshesAt( self, ix, iy, iz ):
        count = len( VoxelGrid.Tetras[0] )
        return [self.TetrahedronMeshAt( ix, iy, iz, index )
            for index in range( count )]

    def InterpolateEdge( self, p, u, q, v, w ):
        t = ( w - u ) / ( v - u )
        return p * ( 1 - t ) + q * t

    Windings = [[0, 0, 0, 0], #  0 -   - 15
                [0, 1, 2, 3], #  1 - T - 14
                [1, 0, 3, 2], #  2 - T - 13
                [0, 1, 2, 3], #  3 - Q - 12

                [2, 3, 0, 1], #  4 - T - 11
                [3, 1, 2, 0], #  5 - Q - 10
                [2, 1, 3, 0], #  6 - Q - 9
                [3, 0, 1, 2], #  7 - T - 8

                [3, 2, 1, 0], #  8 - T - 7
                [1, 2, 3, 0], #  9 - Q - 6
                [0, 2, 1, 3], # 10 - Q - 5
                [2, 1, 0, 3], # 11 - T - 4

                [3, 2, 1, 0], # 12 - Q - 3
                [1, 2, 3, 0], # 13 - T - 2
                [0, 3, 2, 1], # 14 - T - 1
                [0, 0, 0, 0]] # 15 -   - 0

    Polygons = [0, 0, 0, 1,
                0, 1, 1, 0,
                0, 1, 1, 0,
                1, 0, 0, 0]

    def IsoSurfaceAt( self, value ):
        mesh = Mesh( )
        for ix, iy, iz in self.VoxelIndices( ):
            for points, values in self.TetrahedraAt( ix, iy, iz ):
                case = self.Scenario( points, values, value )
                if( case == 0 or case == 15 ): continue
                self.AddFace( mesh, case, points, values, value )
        return mesh

    def Scenario( self, points, values, value ):
        return ( ( 1 if( values[0] < value ) else 0 ) |
                 ( 2 if( values[1] < value ) else 0 ) |
                 ( 4 if( values[2] < value ) else 0 ) |
                 ( 8 if( values[3] < value ) else 0 ) )

    def AddFace( self, mesh, case, points, values, value ):
        points, values = self.Permutate( case, points, values )
        if( VoxelGrid.Polygons[case] ):
            return self.AddQuad( mesh, points, values, value )
        return self.AddTriangle( mesh, points, values, value )

    def Permutate( self, case, points, values ):
        i, j, k, l = VoxelGrid.Windings[case]
        return ( [points[i], points[j], points[k], points[l]],
                 [values[i], values[j], values[k], values[l]] )

    def AddQuad( self, mesh, points, values, value ):
        a, b, c, d = points
        A, B, C, D = values

        e = self.InterpolateEdge( d, D, a, A, value )
        f = self.InterpolateEdge( d, D, b, B, value )
        g = self.InterpolateEdge( c, C, b, B, value )
        h = self.InterpolateEdge( c, C, a, A, value )

        mesh.Faces.AddFace(
            mesh.Vertices.Add( e ),
            mesh.Vertices.Add( f ),
            mesh.Vertices.Add( g ),
            mesh.Vertices.Add( h ) )
        return self

    def AddTriangle( self, mesh, points, values, value ):
        a, b, c, d = points
        A, B, C, D = values

        e = self.InterpolateEdge( a, A, b, B, value )
        f = self.InterpolateEdge( a, A, c, C, value )
        g = self.InterpolateEdge( a, A, d, D, value )

        mesh.Faces.AddFace(
            mesh.Vertices.Add( e ),
            mesh.Vertices.Add( f ),
            mesh.Vertices.Add( g ) )
        return self

    def ImplicitSurface( self, function ):
        for ix, iy, iz in self.PointIndices( ):
            point = self.PointAt( ix, iy, iz )
            self.Values[ix][iy][iz] = function(
                 point.X, point.Y, point.Z )
        return self

    def SignedDistanceField( self, subject ):
        for ix, iy, iz in self.PointIndices( ):
            point = self.PointAt( ix, iy, iz )
            self.Values[ix][iy][iz] = (
                subject.SignedDistance( point ) )

class Sphere:
    def __init__( self, origin, radius ):
        self.Origin = origin
        self.Radius = radius

    def SignedDistance( self, point ):
        o, r = self.Origin, self.Radius
        return o.DistanceTo( point ) - r

class Cylinder:
    def __init__( self, origin, direction, radius ):
        self.Origin    = origin
        self.Direction = direction
        self.Radius    = radius

    def SignedDistance( self, point ):
        o, u, r = self.Origin, self.Direction, self.Radius
        return point.DistanceTo( o + u * ( point - o ) * u ) - r

class BooleanUnion:
    def __init__( self, subjects ):
        self.Subjects = subjects

    def SignedDistance( self, point ):
        distances = [subject.SignedDistance( point )
            for subject in self.Subjects]
        return min( distances )

class BooleanIntersection:
    def __init__( self, subjects ):
        self.Subjects = subjects

    def SignedDistance( self, point ):
        distances = [subject.SignedDistance( point )
            for subject in self.Subjects]
        return max( distances )

class BooleanDifference:
    def __init__( self, subject, _object ):
        self.Subject = subject
        self.Object  = _object

    def SignedDistance( self, point ):
        subject = self.Subject.SignedDistance( point )
        _object = self.Object .SignedDistance( point )
        return max( subject, -_object )

class Offset:
    def __init__( self, object_, distance ):
        self.Object   = object_
        self.Distance = distance

    def SignedDistance( self, point ):
        distance = self.Object.SignedDistance( point )
        return distance - self.Distance

class Shell:
    def __init__( self, object_, thickness ):
        self.Object    = object_
        self.Thickness = thickness

    def SignedDistance( self, point ):
        distance = self.Object.SignedDistance( point )
        return abs( distance ) - self.Thickness

class Blend:
    def __init__( self, subjects, kappa ):
        self.Subjects = subjects
        self.Kappa    = kappa * 4.0
        self.KappaInv = 1.0 / self.Kappa

    def QuadraticMin( self, a, b ):
        k, i = self.Kappa, self.KappaInv
        h = max( k - abs( a - b ), 0.0 ) * i
        return min( a, b ) - h *  h * k * 0.25

    def SignedDistance( self, point ):
        distances = [subject.SignedDistance( point )
            for subject in self.Subjects]
        minimum = distances[0]
        for distance in distances[1:]:
            minimum = self.QuadraticMin(
                minimum, distance )
        return minimum

voxels = VoxelGrid(
    Point3d( -Domain, -Domain, -Domain ),
    Point3d(  Domain,  Domain,  Domain ),
    [Resolution, Resolution, Resolution] )

voxels.SignedDistanceField( Blend( [
    Sphere( Point3d( -0.6, 0.0, 0.0 ), 0.4 ),
    Sphere( Point3d(  0.6, 0.0, 0.0 ), 0.4 ),
    Sphere( Point3d(  0.0, 0.5, 0.5 ), 0.4 )], 0.15 ) )

surface = voxels.IsoSurfaceAt( 0 )
print( surface.Vertices.Count, surface.Faces.Count )

Geometry = [surface]

Show( surface )