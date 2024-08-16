 
        function showClusterLoc(obj,idx)
            assert(~isempty(obj.clustLoc),'Please run getAllClusterLocalization');
            
            locData = obj.clustLoc{idx};
            
            pos = norm([locData.x,locData.y]);
            int = locData.int;
            angle = locData.angle;
            
            figure
            subplot(1,3,1)
            scatter(pos,int,10,'filled')
            axis square
            box on
            title('Intensity vs Position')
            
            subplot(1,3,2)
            scatter(pos,angle,10,'filled')
            axis square
            box on
            title('Position vs angle')
            
            subplot(1,3,3)
            scatter(int,angle,10,'filled')
            axis square
            box on
            title('Intensity vs Angle')
            
        end