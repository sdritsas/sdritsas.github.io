export function print( args ) {
    console.log( args );
}
export class Vec3 {
    constructor({ X = 0.0,
                  Y = 0.0,
                  Z = 0.0 }) {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
    }

    static get O( ) { return new Vec3({ X: 0.0, Y: 0.0, Z: 0.0 }); }
    static get X( ) { return new Vec3({ X: 1.0, Y: 0.0, Z: 0.0 }); }
    static get Y( ) { return new Vec3({ X: 0.0, Y: 1.0, Z: 0.0 }); }
    static get Z( ) { return new Vec3({ X: 0.0, Y: 0.0, Z: 1.0 }); }

    ToStrings( digits = 3, padding = 8 ) {
        const x = this.X.toFixed( digits ).padStart( padding, ' ' );
        const y = this.Y.toFixed( digits ).padStart( padding, ' ' );
        const z = this.Z.toFixed( digits ).padStart( padding, ' ' );
        return [x, y, z];
    }
    ToRowString( ) {
        return this.ToStrings( ).join( '' );
    }
    WithX( x ) {
        this.X = x;
        return this;
    }
    WithY( y ) {
        this.Y = y;
        return this;
    }
    WithZ( z ) {
        this.Z = z;
        return this;
    }
    Distance( that ) {
        const dx = this.X - that.X;
        const dy = this.Y - that.Y;
        const dz = this.Z - that.Z;
        return Math.sqrt(
            dx * dx +
            dy * dy +
            dz * dz );
    }
    DistanceXY( that ) {
        const dx = this.X - that.X;
        const dy = this.Y - that.Y;
        return Math.sqrt(
            dx * dx +
            dy * dy );
    }
    Duplicate( ) {
        return new Vec3({
            X: this.X,
            Y: this.Y,
            Z: this.Z });
    }
    With( that ) {
        this.X = that.X;
        this.Y = that.Y;
        this.Z = that.Z;
        return this;
    }
    static Add( u, v ) {
        return new Vec3({
            X: u.X + v.X,
            Y: u.Y + v.Y,
            Z: u.Z + v.Z });
    }
    Translate( u ) {
        this.X += u.X;
        this.Y += u.Y;
        this.Z += u.Z;
        return this;
    }
    static Sub( u, v ) {
        return new Vec3({
            X: u.X - v.X,
            Y: u.Y - v.Y,
            Z: u.Z - v.Z });
    }
    Scale( s ) {
        this.X *= s;
        this.Y *= s;
        this.Z *= s;
        return this;
    }
    Dot( that ) {
        return this.X * that.X +
               this.Y * that.Y +
               this.Z * that.Z;
    }
    static Mul( u, s ) {
        return new Vec3({
            X: u.X * s,
            Y: u.Y * s,
            Z: u.Z * s });
    }
    static Cross( u, v ) {
        return new Vec3({
            X: u.Y * v.Z - u.Z * v.Y,
            Y: u.Z * v.X - u.X * v.Z,
            Z: u.X * v.Y - u.Y * v.X });
    }
    get Length( ) {
        return Math.sqrt(
            this.X * this.X +
            this.Y * this.Y +
            this.Z * this.Z );
    }
    get Ortho( ) {
        return new Vec3({
            X: this.Y,
            Y:-this.X,
            Z: this.Z });
    }
    Normalize( scale = 1.0 ) {
        let length = this.Length;
        length = length == 0.0 ?
            0.0 : scale / length;
        this.X *= length;
        this.Y *= length;
        this.Z *= length;
        return this;
    }
}
export class Line {
    constructor({ Source = null,
                  Target = null }) {
        this.Source = Source ?? Vec3.O;
        this.Target = Target ?? Vec3.X;
    }
}
export class Vec4 {
    constructor({ X = 0.0,
                  Y = 0.0,
                  Z = 0.0,
                  W = 0.0 }) {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        this.W = W;
    }

    ToList( ) {
        return [this.X, this.Y, this.Z, this.W];
    }

    ToStrings( digits = 3, padding = 8 ) {
        const x = this.X.toFixed( digits ).padStart( padding, ' ' );
        const y = this.Y.toFixed( digits ).padStart( padding, ' ' );
        const z = this.Z.toFixed( digits ).padStart( padding, ' ' );
        const w = this.W.toFixed( digits ).padStart( padding, ' ' );
        return [x, y, z, w]; }
}
export class Mat4 {
    constructor({ xx = 1.0, yx = 0.0, zx = 0.0, wx = 0.0,
                  xy = 0.0, yy = 1.0, zy = 0.0, wy = 0.0,
                  xz = 0.0, yz = 0.0, zz = 0.0, wz = 0.0,
                  xw = 0.0, yw = 0.0, zw = 0.0, ww = 1.0 }) {
        this.X = new Vec4({ X: xx, Y: xy, Z: xz, W: xw });
        this.Y = new Vec4({ X: yx, Y: yy, Z: yz, W: yw });
        this.Z = new Vec4({ X: zx, Y: zy, Z: zz, W: zw });
        this.W = new Vec4({ X: wx, Y: wy, Z: wz, W: ww });
    }
    static get Identity( ) {
        return new Mat4({
            xx: 1.0, yx: 0.0, zx: 0.0, wx: 0.0,
            xy: 0.0, yy: 1.0, zy: 0.0, wy: 0.0,
            xz: 0.0, yz: 0.0, zz: 1.0, wz: 0.0,
            xw: 0.0, yw: 0.0, zw: 0.0, ww: 1.0 });
    }
    static Basis( o, x, y, z ) {
        return new Mat4({
            xx: x.X, yx: y.X, zx: z.X, wx: o.X,
            xy: x.Y, yy: y.Y, zy: z.Y, wy: o.Y,
            xz: x.Z, yz: y.Z, zz: z.Z, wz: o.Z,
            xw: x.W, yw: y.W, zw: z.W, ww: o.W });
    }
    static Translate( tx, ty, tz ) {
        return new Mat4({
            xx: 1.0, yx: 0.0, zx: 0.0, wx: tx,
            xy: 0.0, yy: 1.0, zy: 0.0, wy: ty,
            xz: 0.0, yz: 0.0, zz: 1.0, wz: tz,
            xw: 0.0, yw: 0.0, zw: 0.0, ww: 1.0 });
    }
    static Scale( sx, sy, sz ) {
        return new Mat4({
            xx:  sx, yx: 0.0, zx: 0.0, wx: 0.0,
            xy: 0.0, yy:  sy, zy: 0.0, wy: 0.0,
            xz: 0.0, yz: 0.0, zz:  sz, wz: 0.0,
            xw: 0.0, yw: 0.0, zw: 0.0, ww: 1.0 });
    }
    static RotateX( angle ) {
        const sin = Math.sin( angle );
        const cos = Math.cos( angle );
        return new Mat4({
            xx: 1.0, yx: 0.0, zx: 0.0, wx: 0.0,
            xy: 0.0, yy: cos, zy:-sin, wy: 0.0,
            xz: 0.0, yz: sin, zz: cos, wz: 0.0,
            xw: 0.0, yw: 0.0, zw: 0.0, ww: 1.0 });
    }
    static RotateY( angle ) {
        const sin = Math.sin( angle );
        const cos = Math.cos( angle );
        return new Mat4({
            xx: cos, yx: sin, zx: 0.0, wx: 0.0,
            xy: 0.0, yy: 1.0, zy: 0.0, wy: 0.0,
            xz:-sin, yz: cos, zz: 1.0, wz: 0.0,
            xw: 0.0, yw: 0.0, zw: 0.0, ww: 1.0 });
    }
    static RotateZ( angle ) {
        const sin = Math.sin( angle );
        const cos = Math.cos( angle );
        return new Mat4({
            xx: cos, yx:-sin, zx: 0.0,
            xy: sin, yy: cos, zy: 0.0,
            xz: 0.0, yz: 0.0, zz: 1.0 });
    }
    static AxisAngle( axis, angle ) {
        const sin = Math.sin( angle );
        const cos = Math.cos( angle );
        const soc = 1.0 - cos;
        return new Mat4({
            xx: axis.X * axis.X * soc + cos,
            yx: axis.X * axis.Y * soc - axis.Z * sin,
            zx: axis.X * axis.Z * soc + axis.Y * sin,

            xy: axis.Y * axis.X * soc + axis.Z * sin,
            yy: axis.Y * axis.Y * soc + cos,
            zy: axis.Y * axis.Z * soc - axis.X * sin,

            xz: axis.Z * axis.X * soc - axis.Y * sin,
            yz: axis.Z * axis.Y * soc + axis.X * sin,
            zz: axis.Z * axis.Z * soc + cos });
    }
    static LookAt( e, c, u ) {
        let zx = e.X - c.X;
        let zy = e.Y - c.Y;
        let zz = e.Z - c.Z;

        let yx = u.X;
        let yy = u.Y;
        let yz = u.Z;

        let xx = yy * zz - yz * zy;
        let xy = yz * zx - yx * zz;
        let xz = yx * zy - yy * zx;

        yx = zy * xz - zz * xy;
        yy = zz * xx - zx * xz;
        yz = zx * xy - zy * xx;

        let dx = Math.sqrt( xx * xx + xy * xy + xz * xz );
        let dy = Math.sqrt( yx * yx + yy * yy + yz * yz );
        let dz = Math.sqrt( zx * zx + zy * zy + zz * zz );

        dx = Math.abs( dx <= 1e-9 ) ? 0.0 : 1.0 / dx;
        dy = Math.abs( dy <= 1e-9 ) ? 0.0 : 1.0 / dy;
        dz = Math.abs( dz <= 1e-9 ) ? 0.0 : 1.0 / dz;

        xx *= dx; xy *= dx; xz *= dx;
        yx *= dy; yy *= dy; yz *= dy;
        zx *= dz; zy *= dz; zz *= dz;

        var tx = xx * e.X + xy * e.Y + xz * e.Z;
        var ty = yx * e.X + yy * e.Y + yz * e.Z;
        var tz = zx * e.X + zy * e.Y + zz * e.Z;

        return new Mat4({
            xx: xx, yx: xy, zx: xz, wx:-tx,
            xy: yx, yy: yy, zy: yz, wy:-ty,
            xz: zx, yz: zy, zz: zz, wz:-tz,
            xw:  0, yw: 0,  zw: 0,  ww:  1 });
    }
    static Ortho( xmin = -1, xmax = 1,
                  ymin = -1, ymax = 1,
                  zmin = -1, zmax = 1 ) {
        const dx = 1.0 / ( xmin - xmax );
        const dy = 1.0 / ( ymin - ymax );
        const dz = 1.0 / ( zmin - zmax );

        const xx =-2.0 * dx;
        const yy =-2.0 * dy;
        const zz = 2.0 * dz;

        const tx = ( xmax + xmin ) * dx;
        const ty = ( ymax + ymin ) * dy;
        const tz = ( zmax + zmin ) * dz;

        return new Mat4({
            xx: xx, yx:  0, zx:  0, wx: tx,
            xy:  0, yy: yy, zy:  0, wy: ty,
            xz:  0, yz:  0, zz: zz, wz: tz,
            xw:  0, yw:  0, zw:  0, ww:  1 });
    }
    ToString( ) {
        const [sxx, sxy, sxz, sxw] = this.X.ToStrings( );
        const [syx, syy, syz, syw] = this.Y.ToStrings( );
        const [szx, szy, szz, szw] = this.Z.ToStrings( );
        const [swx, swy, swz, sww] = this.W.ToStrings( );
        return `${sxx} ${syx} ${szx} ${swx}\n` +
               `${sxy} ${syy} ${szy} ${swy}\n` +
               `${sxz} ${syz} ${szz} ${swz}\n` +
               `${sxw} ${syw} ${szw} ${sww}`;
    }
    Multiply( that ) {
        const [sxx, sxy, sxz, sxw] = this.X.ToList( );
        const [syx, syy, syz, syw] = this.Y.ToList( );
        const [szx, szy, szz, szw] = this.Z.ToList( );
        const [swx, swy, swz, sww] = this.W.ToList( );

        const [txx, txy, txz, txw] = that.X.ToList( );
        const [tyx, tyy, tyz, tyw] = that.Y.ToList( );
        const [tzx, tzy, tzz, tzw] = that.Z.ToList( );
        const [twx, twy, twz, tww] = that.W.ToList( );

        const rxx = sxx * txx + syx * txy + szx * txz + swx * txw;
        const ryx = sxx * tyx + syx * tyy + szx * tyz + swx * tyw;
        const rzx = sxx * tzx + syx * tzy + szx * tzz + swx * tzw;
        const rwx = sxx * twx + syx * twy + szx * twz + swx * tww;

        const rxy = sxy * txx + syy * txy + szy * txz + swy * txw;
        const ryy = sxy * tyx + syy * tyy + szy * tyz + swy * tyw;
        const rzy = sxy * tzx + syy * tzy + szy * tzz + swy * tzw;
        const rwy = sxy * twx + syy * twy + szy * twz + swy * tww;

        const rxz = sxz * txx + syz * txy + szz * txz + swz * txw;
        const ryz = sxz * tyx + syz * tyy + szz * tyz + swz * tyw;
        const rzz = sxz * tzx + syz * tzy + szz * tzz + swz * tzw;
        const rwz = sxz * twx + syz * twy + szz * twz + swz * tww;

        const rxw = sxw * txx + syw * txy + szw * txz + sww * txw;
        const ryw = sxw * tyx + syw * tyy + szw * tyz + sww * tyw;
        const rzw = sxw * tzx + syw * tzy + szw * tzz + sww * tzw;
        const rww = sxw * twx + syw * twy + szw * twz + sww * tww;

        return new Mat4({
            xx: rxx, yx: ryx, zx: rzx, wx: rwx,
            xy: rxy, yy: ryy, zy: rzy, wy: rwy,
            xz: rxz, yz: ryz, zz: rzz, wz: rwz,
            xw: rxw, yw: ryw, zw: rzw, ww: rww });
    }
    Point( u ) {
        const [xx, yx, zx] = [this.X.X, this.Y.X, this.Z.X];
        const [xy, yy, zy] = [this.X.Y, this.Y.Y, this.Z.Y];
        const [xz, yz, zz] = [this.X.Z, this.Y.Z, this.Z.Z];
        const [ox, oy, oz] = [this.W.X, this.W.Y, this.W.Z];
        return new Vec3({
            X: ox + xx * u.X + yx * u.Y + zx * u.Z,
            Y: oy + xy * u.X + yy * u.Y + zy * u.Z,
            Z: oz + xz * u.X + yz * u.Y + zz * u.Z });
    }
    Vector( u ) {
        const [xx, yx, zx] = [this.X.X, this.Y.X, this.Z.X];
        const [xy, yy, zy] = [this.X.Y, this.Y.Y, this.Z.Y];
        const [xz, yz, zz] = [this.X.Z, this.Y.Z, this.Z.Z];
        return new Vec3({
            X: xx * u.X + yx * u.Y + zx * u.Z,
            Y: xy * u.X + yy * u.Y + zy * u.Z,
            Z: xz * u.X + yz * u.Y + zz * u.Z });
    }
}
export class Palette {
    static Dark      = 0;
    static Light     = 1;
    static Scheme    = Palette.Dark;

    static Foreground = [ ' #ffffff', ' #000000' ];
    static Background = [ ' #000000', ' #ffffff' ];

    static AxisX      = [ ' #ff4500', ' #873213' ];
    static AxisY      = [ ' #00ff45', ' #0b9e32' ];
    static AxisZ      = [ ' #0045ff', ' #0333b6' ];

    static NodeEdgeHigh = [ ' #ffcc00', ' #c08300' ];
    static NodeFillHigh = [ ' #ffe682', ' #c08300' ];
    static NodeEdgeBase = [ ' #888888', ' #888888' ];
    static NodeFaceBase = [ ' #BBBBBB', ' #BBBBBB' ];
}
export class FillStyle {
    WithColor( color ) {
        this.Color = color;
        return this;
    }
    static get Default( ) {
        return new FillStyle( )
            .WithColor( Palette.Background );
    }
}
export class StrokeSymbol {
    static Plain  = 'plain';
    static Arrow  = 'arrow';
    static Limit  = 'limit';
    WithSymbol( symbol = StrokeSymbol.Plain ) {
        this.Symbol = symbol;
        return this;
    }
    WithHandle( handle = false ) {
        this.Handle = handle;
        return this;
    }
    get IsPlain( ) {
        return this.Symbol == StrokeSymbol.Plain;
    }
    get IsArrow( ) {
        return this.Symbol == StrokeSymbol.Arrow;
    }
    get IsLimit( ) {
        return this.Symbol == StrokeSymbol.Limit;
    }
    static get Default( ) {
        return new StrokeSymbol( )
            .WithSymbol( )
            .WithHandle( );
    }
    static get Vector( ) {
        return StrokeSymbol.Default
            .WithSymbol( StrokeSymbol.Arrow );
    }
    static get VectorHandle( ) {
        return StrokeSymbol.Vector
            .WithHandle( true );
    }
}
export class StrokePattern {
    static None = [];
    static Dash = [5,5]
}
export class StrokeWidth {
    static HairLine = 1;
    static BaseLine = 2;
    static Vec3Line = 5;
}
export class StrokeStyle {
    WithColor( color = Palette.Foreground ) {
        this.Color = color;
        return this;
    }
    WithWidth( width = 1.0 ) {
        this.Width = width;
        return this;
    }
    WithPattern( pattern = StrokePattern.None ) {
        this.Pattern = pattern;
        return this;
    }
    WithSource( source = null )
    {
        this.Source = source ?? StrokeSymbol.Default;
        return this;
    }
    WithTarget( target = null )
    {
        this.Target = target ?? StrokeSymbol.Default;
        return this;
    }
    static get Default( ) {
        return new StrokeStyle( )
            .WithColor( )
            .WithWidth( )
            .WithPattern( )
            .WithSource( )
            .WithTarget( );
    }
}
export class Layers {
    static Foreground = 1000;
    static Background = 0;
    static HiddenLine =-1000;
}
export class Style {
    WithVisible( visible = true ) {
        this.Visible = visible;
        return this;
    }
    WithLayer( layer = Layers.Foreground ) {
        this.Layer = layer;
        return this;
    }
    WithInteractive( interactive = false ) {
        this.Interactive = interactive;
        if( interactive ) {
            this.WithHighlightedStroke( )
                .WithHighlightedFill( );
        }
        this.Highlighted = false;
        return this;
    }
    WithStroke( stroke = null ) {
        this.RegularStroke = stroke ?? StrokeStyle.Default;
        return this;
    }
    WithHighlightedStroke( stroke = null ) {
        this.HighlightedStroke = stroke ?? StrokeStyle.Default;
        return this;
    }
    get Stroke( ) {
        return this.Highlighted ?
            this.HighlightedStroke :
            this.RegularStroke;
    }
    WithFill( fill = null ) {
        this.RegularFill = fill ?? FillStyle.Default;
        return this;
    }
    WithHighlightedFill( fill = null ) {
        this.HighlightedFill = fill ?? FillStyle.Default;
        return this;
    }
    get Fill( ) {
        return this.Highlighted ?
            this.HighlightedFill :
            this.RegularFill;
    }
    WithDefaults( ) {
        return this.WithVisible( )
                   .WithLayer( )
                   .WithInteractive( )
                   .WithStroke( )
                   .WithFill( );
    }
    static get Default( ) {
        return new Style( )
            .WithDefaults( );
    }
}
export class NodeStyle extends Style {
    static DefaultRadius = 5;
    static DefaultWidth  = 4;
    WithRadius( radius = NodeStyle.DefaultRadius ) {
        this.Radius = radius;
        return this;
    }
    static get Default( ) {
        return new NodeStyle( )
            .WithDefaults( )
            .WithFill( FillStyle.Default
                .WithColor( Palette.NodeFaceBase ) )
            .WithStroke( StrokeStyle.Default
                .WithColor( Palette.NodeEdgeBase )
                .WithWidth( NodeStyle.DefaultWidth ) )
            .WithRadius( );
    }
    static get Hidden( ) {
        return NodeStyle.Default
            .WithVisible( false );
    }
    static get Interactive( ) {
        return NodeStyle.Default
            .WithInteractive( true )
            .WithHighlightedFill( FillStyle.Default
                .WithColor( Palette.NodeFillHigh ) )
            .WithHighlightedStroke( StrokeStyle.Default
                .WithColor( Palette.NodeEdgeHigh )
                .WithWidth( NodeStyle.DefaultWidth ) );
    }
}
export class EdgeStyle extends Style {
    static get Default( ) {
        return new EdgeStyle( )
            .WithDefaults( );
    }
    static Line( color,
                 width = StrokeWidth.BaseLine,
                 dash  = StrokePattern.None ) {
        return new EdgeStyle( )
            .WithDefaults( )
            .WithLayer( Layers.Background )
            .WithFill( FillStyle.Default
                .WithColor( color ) )
            .WithStroke( StrokeStyle.Default
                .WithColor( color )
                .WithPattern( dash )
                .WithWidth( width ) );
    }
    static Dash( color,
                 width = StrokeWidth.BaseLine ) {
        return EdgeStyle.Line( color, width, StrokePattern.Dash )
            .WithLayer( Layers.HiddenLine );
    }
    static get DashX( ) {
        return EdgeStyle.Dash( Palette.AxisX )
    }
    static get DashY( ) {
        return EdgeStyle.Dash( Palette.AxisY )
    }
    static get DashZ( ) {
        return EdgeStyle.Dash( Palette.AxisZ )
    }
    static Vector( color, handle = false ) {
        return new EdgeStyle( )
            .WithDefaults( )
            .WithLayer( Layers.Background )
            .WithFill( FillStyle.Default
                .WithColor( color ) )
            .WithStroke( StrokeStyle.Default
                .WithColor( color )
                .WithTarget( handle ?
                    StrokeSymbol.VectorHandle :
                    StrokeSymbol.Vector )
                .WithWidth( StrokeWidth.Vec3Line ) );
    }
    static get AxisX( ) {
        return EdgeStyle.Vector(
            Palette.AxisX, false );
    }
    static get AxisY( ) {
        return EdgeStyle.Vector(
            Palette.AxisY, false );
    }
    static get AxisZ( ) {
        return EdgeStyle.Vector(
            Palette.AxisZ, false );
    }
    static get HandleAxisX( ) {
        return EdgeStyle.Vector(
            Palette.AxisX, true );
    }
    static get HandleAxisY( ) {
        return EdgeStyle.Vector(
            Palette.AxisY, true );
    }
    static get HandleAxisZ( ) {
        return EdgeStyle.Vector(
            Palette.AxisZ, true );
    }
}
export class TextStyle extends Style { 
    static get Default( ) {
        return new TextStyle( )
            .WithDefaults( );
    }
    WithTextFont( font = 'times' ) {
        this.TextFont = font;
        return this;
    }
    WithTextSize( size = 10 ) {
        this.TextSize = size;
        return this;
    }
    WithTextColor( color = '#ffffff' ) {
        this.TextColor = color;
        return this;
    }
}
export const Reactive = ( Type ) => class extends Type {
    WithDefaults( ) {
        return this.WithName( )
                   .WithStyle( )
                   .WithDependencies( )
                   .WithUpdate( );
    }
    WithDependencies( dependencies = null ) {
        this.Dependencies = dependencies ?? [];
        return this;
    }
    WithUpdate( update = null ) {
        this.Update = update ?? ( ( ) => { } );
        return this;
    }
    WithName( name = '' ) {
        this.Name = name;
        return this;
    }
    WithStyle( style = null ) {
        this.Style = style;
        return this;
    }
    Change( that ) {
        this.With( that );
        this.Notify( );
    }
    Notify( ) {
        for( const dependency of this.Dependencies ) {
            dependency.Update( );
        }
    }
}
export function Computed ({ Class  = null,
                            Method = null,
                            Inputs = null }) {
    const result = Method( ...Inputs );
    const self = Class.From( result );
    self.Properties( ...arguments );
    self.Update = ( ) => {
        self.Change( Method( ...Inputs ) );
    };
    for( const input of Inputs ) {
        input.Deps.push( self );
    }
    return self;
}
export class Node extends Reactive( Vec3 ) {
    constructor( x = 0.0,
                 y = 0.0,
                 z = 0.0 ) {
        super({ X: x, Y: y, Z: z });
        this.WithDefaults( )
            .WithStyle( NodeStyle.Default );
    }
    static FromVec3( u ) {
        return new Node( u.X, u.Y, u.Z );
    }
    ToScreen( viewport ) {
        return viewport.Transform( this );
    }
    Overlaps( viewport, mouse ) {
        const point = this.ToScreen( viewport );
        return point.DistanceXY( mouse ) <= this.Style.Radius;
    }
    Paint( viewport ) {
        viewport.Style( this.Style );
        viewport.PaintPath( );
        viewport.PaintCircle(
            this.ToScreen( viewport ),
            this.Style.Radius );
        viewport.PaintBoth( );
    }
}
export class Edge extends Reactive( Line ) {
    constructor( source, target  ) {
        super({ Source: source, Target: target });
        this.WithDefaults( )
            .WithStyle( EdgeStyle.Default );
    }
    PaintArrow( view, o, x ) {
        const y = x.Ortho;
        const points = [o];
        for( const [sx, sy] of [[2, 0], [2, 1],
                                [1, 1],
                                [1, 2], [0, 2]] ) {
            points.push( Vec3.Add( o,
                Vec3.Add( Vec3.Mul( x, sx ),
                          Vec3.Mul( y, sy ) ) ) );
        }
        view.PaintPath( );
        view.PaintPoly( points );
        view.PaintFill( );
    }
    Paint( viewport ) {
        const source = viewport.Transform( this.Source );
        const target = viewport.Transform( this.Target );

        const distance = source.Distance( target );
        if( distance <= 1 ) return;

        const u = Vec3.Sub( target, source )
                      .Scale( 1.0 / distance );
        const v = u.Ortho;

        if( this.Style.Stroke.Source.Handle ) {
            //const radius = this.Source.Style.Radius;
            source.Translate( Vec3.Mul( u,
                NodeStyle.DefaultRadius ) );
        }

        if( this.Style.Stroke.Target.Handle )
            target.Translate( Vec3.Mul( u,
                -NodeStyle.DefaultRadius ) );

        viewport.Style( this.Style );

        if( distance > 15 && this.Style.Stroke.Target.IsArrow ) {
            this.PaintArrow( viewport, target,
                Vec3.Sub( v, u ).Normalize( 10 ) );
            target.Translate( Vec3.Mul( u, -10 ) );
        }

        viewport.PaintPath( );
        viewport.PaintLine( source, target );
        viewport.PaintEdge( );
    }
}
export class Camera {
    constructor( source, target, normal, scale ) {
        this.Source = source;
        this.Target = target;
        this.Normal = normal;
        this.Scale  = scale;
    }
    static get Default( ) {
        return new Camera(
            new Vec3({ X: 0.0, Y: 0.0, Z: 1.0 }),
            new Vec3({ X: 0.0, Y: 0.0, Z: 0.0 }),
            new Vec3({ X: 0.0, Y: 1.0, Z: 0.0 }), 10.0 );
    }
    get Frame( ) {
        const o = this.Source;
        const z = Vec3.Sub( this.Target, this.Source ).Normalize( );
        const y = this.Normal;
        const x = Vec3.Cross( y, z );
        return [o, x, y, z]
    }
    get Direction( ) {
        return Vec3.Sub( this.Target,
                         this.Source );
    }
    get Distance( ) {
        return this.Source.Distance( this.Target );
    }
    get Matrix( ) {
        return Mat4.LookAt( this.Source,
                            this.Target,
                            this.Normal );
    }
}
export class Model {
    constructor( ) {
        this.Elements = [];
    }
    Sort( ) {
        this.Elements.sort( ( a, b ) =>
            a.Style.Layer - b.Style.Layer );
        return this;
    }
    Add( item ) {
        this.Elements.push( item );
        return item;
    }
    AddNode( x, y, z ) {
        return this.Add(
            new Node( x, y, z ) );
    }
    AddEdge( source, target ) {
        return this.Add(
            new Edge( source, target ) );
    }
    AddLive( inputs, method ) {
        const result = method( ...inputs );
        if( result instanceof Vec3 ) {
            const node = Node.FromVec3( result );
            node.Update = ( ) => {
                node.Change( method( ...inputs ) );
            };
            for( const input of inputs ) {
                input.Dependencies.push( node );
            }
            return this.Add( node );
        }
        return null;
    }
    // AddWorld( ) {
    //     this.AddArrow( Vec3.O, Vec3.X, { Edge: System.ColorAxisX   } );
    //     this.AddArrow( Vec3.O, Vec3.Y, { Edge: System.ColorAxisY } );
    //     this.AddArrow( Vec3.O, Vec3.Z, { Edge: System.ColorAxisZ  } );
    //     this.AddPoint( Vec3.O );
    // }
}
export class Viewport {
    constructor({ Model  = null,
                  Figure = 'none',
                  Rotate = true,
                  Width  = 0,
                  Height = 256,
                  World  = true }) {
        this.Model  = Model ?? new Model( );
        this.Model.Sort( );
        this.Figure = Figure;
        this.Rotate = Rotate;
        this.Width  = Width;
        this.Height = Height;
        this.Matrix = Mat4.Identity;
        this.Camera = Camera.Default;
        if( World ) this.AddWorld( );

        this.Script = document.getElementById( this.Figure );
        this.Width = this.Script.parentElement.clientWidth;

        this.Canvas = document.createElement( 'canvas' );
        this.Canvas.setAttribute( 'id', `${this.Figure}-canvas` );
        this.Canvas.setAttribute( 'style', `width:100%;height:${this.Height}px` );
        this.Canvas.setAttribute( 'width',  this.Width  );
        this.Canvas.setAttribute( 'height', this.Height );
        this.Canvas.setAttribute( 'viewport', this );

        this.MouseStatus = 'none';
        this.MouseCamera = Camera.Default;
        this.MouseSource = new Vec3({ X: 0, Y: 0 });
        this.MouseTarget = new Vec3({ X: 0, Y: 0 });
        this.Canvas.addEventListener( 'mousedown', event => this.MouseDown    ( event ) );
        this.Canvas.addEventListener( 'mouseup',   event => this.MouseUp      ( event ) );
        this.Canvas.addEventListener( 'mousemove', event => this.MouseMove    ( event ) );

        this.Script.parentNode.insertBefore( this.Canvas, this.Script.nextSibling );

        this.Context = this.Canvas.getContext( '2d' );
        this.Context.lineCap  = 'round';
        this.Context.lineJoin = 'round';

        System.Viewports.push( this );
        this.Ready = true;
        return this;
    }
    Clear( ) {
        this.Context.clearRect( 0, 0,
            this.Width, this.Height );
    }
    UpdateMatrix( ) {
        const tx = 0.5 * this.Width;
        const ty = 0.5 * this.Height;

        let dx = tx / ty;
        let dy = 1.0;
        if( dx < 1.0 ) {
            dy = 1.0 / dx;
            dx = 1.0; }
        dx = dx / tx;
        dy =-dy / ty;

        const P = Mat4.Ortho( -dx, dx, -dy, dy );
        const T = Mat4.Translate( tx, ty, 0 );
        const C = this.Camera.Matrix;
        this.Matrix = T.Multiply( P.Multiply( C ) );
        //console.log( 'Matrix\n' + this.Matrix.ToString( ) );
    }
    Transform( point ) {
        return this.Matrix.Point( point ).WithZ( 0.0 );
    }
    PaintPath( ) {
        this.Context.beginPath( );
    }
    PaintEdge( ) {
        this.Context.stroke( );
    }
    PaintFill( ) {
        this.Context.fill( );
    }
    PaintBoth( ) {
        this.Context.stroke( );
        this.Context.fill( );
    }
    PaintMove( point ) {
        this.Context.moveTo( point.X, point.Y );
    }
    PaintLine( source, target ) {
        this.Context.moveTo( source.X, source.Y );
        this.Context.lineTo( target.X, target.Y );
    }
    PaintPoly( points ) {
        let point = points[0];
        this.Context.moveTo( point.X, point.Y );
        for( let index = 1; index < points.length; index++ ) {
            point = points[index];
            this.Context.lineTo( point.X, point.Y );
        }
    }
    PaintCircle( origin, radius ) {
        this.Context.ellipse( origin.X, origin.Y,
            radius, radius, 0.0, 0.0, Math.PI * 2.0 );
        return this;
    }
    Style( style ) {
        this.Context.fillStyle   = style.Fill.Color[Palette.Scheme];
        this.Context.strokeStyle = style.Stroke.Color[Palette.Scheme];
        this.Context.lineWidth   = style.Stroke.Width;
        this.Context.setLineDash( style.Stroke.Pattern );
        return this;
    }
    Paint( ) {
        this.UpdateMatrix( );
        for( const element of this.Model.Elements ) {
            if( !element.Style.Visible ) continue;
            element.Paint( this );
        }
    }
    Repaint( ) {
        this.Clear( );
        this.Paint( );
    }
    MouseDownNone( ) {
        for( const element of this.Model.Elements ) {
            if( element.Style.Interactive ) {
                if( element.Overlaps( this, this.MouseSource ) ) {
                    this.MouseObject = [element,
                                        element.Duplicate( )];
                    this.MouseStatus = 'drag';
                    return false;
                }
            }
        }
        if( this.Rotate ) {
            this.MouseCamera.Source.With( this.Camera.Source );
            this.MouseCamera.Normal.With( this.Camera.Normal );
            this.MouseStatus = 'turn';
        }
        return false;
    }
    MouseMoveNone( ) {
        let changed = false;
        for( const element of this.Model.Elements ) {
            if( element.Style.Interactive ) {
                if( element.Overlaps( this, this.MouseTarget ) ) {
                    if( !element.Style.Highlighted ) {
                        element.Style.Highlighted = true;
                        changed = true;
                    }
                }
                else
                {
                    if( element.Style.Highlighted ) {
                        element.Style.Highlighted = false;
                        changed = true;
                    }
                }
            }
        }
        if( changed )
            this.Repaint( );
        return false;
    }
    MouseMoveDrag( ) {
        const delta = Vec3.Sub(
            this.MouseTarget,
            this.MouseSource );
        if( delta.Length <= 1 )
            return false;

        const [ o, x, y, z ] = this.Camera.Frame;
        const distance = 2.0 * this.Camera.Distance;

        const sx = distance * this.Width / this.Height;
        const sy = distance * 1.0;

        delta.X *=-sx / this.Width;
        delta.Y *=-sy / this.Height;

        const [element, point] = this.MouseObject;
        element.Change( Vec3.Add( point, Vec3.Add(
            Vec3.Mul( x, delta.X ),
            Vec3.Mul( y, delta.Y ) ) ) );
        this.Repaint( );
    }
    MouseDown( event ) {
        if( !this.Ready ) return true;
        event.preventDefault( );
        event.stopPropagation( );

        this.MouseSource.X = event.offsetX;
        this.MouseSource.Y = event.offsetY;
        switch( this.MouseStatus ) {
            case 'none':
                return this.MouseDownNone( );
            default:
                return false;
        }
    }
    MouseMove( event ) {
        if( !this.Ready ) return true;
        event.preventDefault( );
        event.stopPropagation( );

        this.MouseTarget.X = event.offsetX;
        this.MouseTarget.Y = event.offsetY;
        switch( this.MouseStatus ) {
            case 'none':
                return this.MouseMoveNone( );

            case 'turn': {
                const delta = Vec3.Sub( this.MouseTarget, this.MouseSource );
                if( delta.Length <= 2 ) return false;
                const n = Vec3.Sub( this.MouseCamera.Source, this.Camera.Target );
                let u = Vec3.Cross( n, Vec3.Z );
                if( u.Length <= 1e-3 )
                    u = new Vec3({ X:-1, Y:0, Z:0 });
                const distance = n.Length;
                u.Normalize( );
                n.Normalize( );
                const ax = Math.PI * delta.X / this.Width;
                const ay = Math.PI * delta.Y / this.Height;

                const rotate = Mat4.RotateZ( ax )
                    .Multiply( Mat4.AxisAngle( u, ay ) );

                this.Camera.Source = Vec3.Add( this.Camera.Target,
                    rotate.Vector( n ).Normalize( distance ) );

                this.Camera.Normal = rotate.Vector(
                    this.MouseCamera.Normal ).Normalize( );

                this.Repaint( );
                return false;
            }
            case 'drag':
                return this.MouseMoveDrag( );
            default:
                return false;
        }
    }
    MouseUp( event ) {
        if( !this.Ready ) return true;
        event.preventDefault( );
        event.stopPropagation( );
        switch( this.MouseStatus ) {
            case 'none':
                return false;
            case 'turn':
                this.MouseStatus = 'none';
                return false;
            default:
                this.MouseStatus = 'none';
                return false;
        }
    }
}
export class System {
    static Viewports   = [];

    static OnLoadEvents = []
    static OnLoadEvent( callback = null ) {
        //-- System On Load Events
        //--
        if( callback === null ) {
            window.addEventListener( 'load', event => {
                for( const callback of System.OnLoadEvents ) {
                    callback( event );
                }
            } );
            //-- Color Scheme Event
            //--
            callback = event => {
                for( const element of document.getElementsByClassName( 'md-header__button' ) )
                {
                    if( !element.hasAttribute( 'title' )  ||
                        !element.getAttribute( 'title' )
                                .includes( 'dark' ) ) continue;
                    element.parentElement.addEventListener( 'change', event => {
                            Palette.Scheme = element.hasAttribute( 'hidden' ) ?
                                Palette.Dark : Palette.Light;
                            for( const viewport of System.Viewports )
                                viewport.Repaint( );
                            } );
                    break;
                }
            };
        }
        this.OnLoadEvents.push( callback );
    }
}
System.OnLoadEvent( );