import { System, Model, Viewport,
         NodeStyle, EdgeStyle, Vec3 } from "./viewport.js";
System.OnLoadEvent( ( ) => {
    const model = new Model( );

    const o = model.AddNode( 0.00, 0.00, 0.00 )
                   .WithStyle( NodeStyle.Hidden )
                   .WithName( 'o' );

    const p = model.AddNode( 0.25, 0.75, 0.00 )
                   .WithStyle( NodeStyle.Interactive )
                   .WithName( 'p' );

    const q = model.AddNode( 0.75, 0.25, 0.00 )
                   .WithStyle( NodeStyle.Interactive )
                   .WithName( 'q' );

    const x = model.AddLive( [o, p, q], ( o, p, q ) => {
                        const u = Vec3.Sub( p, o );
                        const v = Vec3.Sub( q, o );
                        return Vec3.Add( o, Vec3.Add( u, v ) ); })
                   .WithStyle( NodeStyle.Hidden )
                   .WithName( 'x' );

    model.AddEdge( o, p )
         .WithStyle( EdgeStyle.HandleAxisX )
         .WithName( 'u' );

    model.AddEdge( o, q )
         .WithStyle( EdgeStyle.HandleAxisY )
         .WithName( 'v' );

    model.AddEdge( o, x )
         .WithStyle( EdgeStyle.AxisZ )
         .WithName( 'w' );

    model.AddEdge( q, x )
         .WithStyle( EdgeStyle.DashX )
         .WithName( 'qx' );

    model.AddEdge( p, x )
         .WithStyle( EdgeStyle.DashY )
         .WithName( 'px' );

    const view = new Viewport({
        Model:  model,
        Figure: 'vector-add',
        Height: 256,
        World:  false,
        Rotate: false } );
    view.Paint( );
} );
